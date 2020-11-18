/*
 * cuts.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"

extern configType config;

int resolveInfeasibility(probType **prob, cellType *cell, bool *newOmegaFlag, int omegaIdx);
int formFeasCut(probType *prob, cellType *cell);
int updtFeasCutPool(numType *num, coordType *coord, cellType *cell);
int checkFeasCutPool(cellType *cell, int lenX);
int addCut2Pool(cellType *cell, oneCut *cut, int lenX, double lb, bool feasCut);

int formSDCut(probType **prob, cellType *cell, dVector Xvect, double lb) {
	oneCut 	*cut;
	dVector piCbarX;
	iVector istar;
	int    	cutIdx;
	bool 	solveSP;

	istar = (iVector) arr_alloc(cell->omega->cnt, int);
	/* Pre-compute pi x bBar - pi x Cbar x x as it is independent of observations */
	piCbarX= arr_alloc(cell->sigma->cnt, double);
	for (int cnt = 0; cnt < cell->sigma->cnt; cnt++)
		piCbarX[cnt] = cell->sigma->vals[cnt].pib - vXv(cell->sigma->vals[cnt].piC, Xvect, prob[1]->coord->CCols, prob[1]->num->cntCcols);

	for ( int obs = 0; obs < cell->omega->cnt; obs++ ) {
		solveSP = false;

		/* Check to see if the observation belongs to the current branch. */
		int cnt = 0;
		while ( cnt < cell->sample->cnt ) {
			if ( obs == cell->sample->omegaIdx[cnt] )
				break;
			cnt++;
		}

		/* A. Determine whether a subproblem will be solved for the observation or the argmax procedure used. */
		switch ( config.SP_SAMPLING ) {
		case Full:
		{
			solveSP = true; break;
		}
		case Binomial:
		{
			double prob;
			prob = scalit(0,1, &config.RUN_SEED[0]);
			if ( prob <= config.SP_FRACTION ) {
				solveSP = true;
			}
			break;
		}
		default:
		{
			if ( cnt < cell->sample->cnt ) {
				solveSP = true;
			}
		}
		}

		if ( solveSP ) {
			/* Option 1: Solve a subproblem */
			/* (a) Construct the subproblem with input observation and master solution, solve the subproblem, and complete stochastic updates */
			if ( (istar[obs] = solveSubprob(prob[1], cell->subprob, Xvect, cell->basis, cell->lambda, cell->sigma, cell->delta,
					config.MAX_ITER, cell->omega, obs, &cell->sample->newOmegaFlag[cnt], cell->k, config.TOLERANCE,
					&cell->spFeasFlag, &cell->sample->newBasisFlag[cnt], &cell->time.subprobIter, &cell->time.argmaxIter) < 0) ){
				errMsg("algorithm", "formSDCut", "failed to solve the subproblem", 0);
				return -1;
			}

			/* (b) Increment the number of subproblems solved during algorithm */
			cell->LPcnt++;

			if ( ! cell->spFeasFlag ) {
				/* (c1) Subproblem is infeasible, resolve infeasibility */
				if ( resolveInfeasibility(prob, cell, &cell->sample->newOmegaFlag[obs], cell->sample->omegaIdx[obs]) ) {
					errMsg("algorithm", "formSDCut", "failed to resolve infeasibility", 0);
					return -1;
				}
			}
			else if ( cell->fcutsPool->cnt > 0 && cell->sample->newOmegaFlag[obs] ) {
				/* (c2) Subproblem is feasible, however a new observation or sigma has been encountered. Therefore, update the feasibility cut
				 * pool and check to see if new feasibility cuts need to be added. */
				if ( formFeasCut(prob[1], cell) ) {
					errMsg("algorithm", "formSDCut", "failed to add new feasibility cuts", 0);
					return -1;
				}
			}
		}
		else {
			/* Option 2: Use the argmax procedure */
			istar[obs] = computeIstar(prob[1]->num, prob[1]->coord, cell->basis, cell->sigma, cell->delta, piCbarX, Xvect,
					cell->omega->vals[obs], obs, 0, cell->basis->cnt);
		}
	}

	/* B. Create an affine lower bound */
	clock_t tic = clock();
	cut = SDCut(prob[1]->num, prob[1]->coord, cell->omega, cell->basis, cell->sigma, cell->delta, istar);
	cell->time.argmaxIter += ((double) (clock()-tic))/CLOCKS_PER_SEC;

#if defined(STOCH_CHECK)
	bool newOmegaFlag = false;
	printf("------\n");
	/* Solve the subproblem to verify if the argmax operation yields a lower bound */
	for ( int cnt = 0; cnt < cell->omega->cnt; cnt++ ) {
		/* (a) Construct the subproblem with input observation and master solution, solve the subproblem, and complete stochastic updates */
		if ( solveSubprob(prob[1], cell->subprob, Xvect, cell->basis, cell->lambda, cell->sigma, cell->delta, config.MAX_ITER,
				cell->omega, cnt, &newOmegaFlag, cell->k, config.TOLERANCE, &cell->spFeasFlag, NULL,
				&cell->time.subprobIter, &cell->time.argmaxIter) < 0 ) {
			errMsg("algorithm", "formSDCut", "failed to solve the subproblem", 0);
			return -1;
		}
	}
	printf("------\n");
#endif

	/* C. Add cut to the structure and master problem  */
	if ( (cutIdx = addCut2Pool(cell, cut, prob[0]->num->cols, lb, false)) < 0) {
		errMsg("algorithm", "formSDCut", "failed to add the new cut to cutsType structure", 0);
		return -1;
	}
	if ( addCut2Master(cell->master, cut, cell->incumbX, prob[0]->num->cols) ) {
		errMsg("algorithm", "formSDCut", "failed to add the new cut to master problem", 0);
		return -1;
	}

	mem_free(istar);
	return cutIdx;
}//END formCut()

oneCut *SDCut(numType *num, coordType *coord, omegaType *omega, basisType *basis, sigmaType *sigma, deltaType *delta, iVector istar) {
	oneCut *cut;
	double alpha, multiplier;
	dVector beta;

	/* allocate memory to hold a new cut */
	cut = newCut(num->prevCols, omega->cnt, omega->sampleSize);

	/* Temporary variables used to update the coefficients */
	alpha = 0.0;
	beta  = (dVector) arr_alloc(num->prevCols + 1, double);

	for (int obs = 0; obs < omega->cnt; obs++) {
		cut->iStar[obs] = istar[obs];

		if ( num->rvdOmCnt > 0 ) {
			for ( int idx = 0; idx <= basis->vals[istar[obs]]->phiLength; idx++ ) {
				int sigmaIdx = basis->vals[istar[obs]]->sigmaIdx[idx];
				int lambdaIdx = sigma->lambdaIdx[sigmaIdx];
				if ( idx == 0 )
					multiplier = 1.0;
				else
					multiplier = omega->vals[obs][coord->rvOffset[2] + basis->vals[istar[obs]]->omegaIdx[idx]];

				/* Start with (Pi x bBar) + (Pi x bomega) + (Pi x Cbar) x X */
				alpha += omega->weights[obs] * multiplier * (sigma->vals[sigmaIdx].pib + delta->vals[lambdaIdx][obs].pib);

				for (int c = 1; c <= num->cntCcols; c++)
					beta[coord->CCols[c]] += omega->weights[obs] * multiplier * sigma->vals[sigmaIdx].piC[c];
				for (int c = 1; c <= num->rvCOmCnt; c++)
					beta[coord->rvCOmCols[c]] += omega->weights[obs] * multiplier * delta->vals[lambdaIdx][obs].piC[c];
			}
		}
		else {
			alpha += sigma->vals[istar[obs]].pib * omega->weights[obs];
			alpha += delta->vals[sigma->lambdaIdx[istar[obs]]][obs].pib * omega->weights[obs];

			for (int c = 1; c <= num->cntCcols; c++)
				beta[coord->CCols[c]] += sigma->vals[istar[obs]].piC[c] * omega->weights[obs];
			for (int c = 1; c <= num->rvCOmCnt; c++)
				beta[coord->rvCols[c]] += delta->vals[sigma->lambdaIdx[istar[obs]]][obs].piC[c] * omega->weights[obs];
		}
	}

	/* Scale the coefficients based on the number of observations */
	cut->alpha = alpha / omega->sampleSize;
	for (int c = 1; c <= num->prevCols; c++)
		cut->beta[c] = beta[c] / omega->sampleSize;
	cut->beta[0] = 1.0;								/* coefficient of eta coloumn */

	mem_free(beta);

	return cut;
}//END SDCut

/* This function loops through a set of cuts and find the highest cut height at the specified position x */
double maxCutHeight(cutsType *cuts, int currSampleSize, dVector xk, int betaLen, double lb) {
	double Sm = -INF, ht = 0.0;
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		ht = cutHeight(cuts->vals[cnt], currSampleSize, xk, betaLen, lb);
		if (Sm < ht) {
			Sm = ht;
		}
	}

	return Sm;
}//END maxCutHeight

/* This function calculates and returns the height of a given cut at a given X.  It includes the k/(k-1) update, but does not include
 * the coefficients due to the cell. */
double cutHeight(oneCut *cut, int currSampleSize, dVector xk, int betaLen, double lb) {
	double height;
	double t_over_k = ((double) cut->numSamples / (double) currSampleSize);

	/* A cut is calculated as alpha - beta x X */
	height = cut->alpha - vXv(cut->beta, xk, NULL, betaLen);

	/* Weight cut based on number of observations used to form it */
	height *= t_over_k;

	/* Updated for optimality cut height*/
	height += (1 - t_over_k) * lb;

	return height;
}//END cutHeight()

/* This function allocates memory for the arrays inside a single cut, and initializes its values accordingly.  The cut structure
 * itself is assumed to be already allocated.  Note, each beta dVector contains room for its one-norm, thought it just gets filled
 * with zero anyway. */
oneCut *newCut(int numX, int numIstar, int numSamples) {
	oneCut *cut;

	cut = (oneCut *) mem_malloc (sizeof(oneCut));
	cut->numSamples = numSamples;
	cut->omegaCnt = numIstar;
	cut->isIncumb = false; 								/* new cut is by default not an incumbent */
	cut->alphaIncumb = 0.0;
	cut->rowNum = -1;

	if ( numIstar > 0 ) {
		if (!(cut->iStar = arr_alloc(numIstar, int)))
			errMsg("allocation", "new_cut", "iStar", 0);
	}
	else
		cut->iStar = NULL;

	if (!(cut->beta = arr_alloc(numX + 1, double)))
		errMsg("allocation", "new_cut", "beta", 0);

	cut->alpha = 0.0;

	cut->name = (cString) arr_alloc(NAMESIZE, char);

	return cut;
}//END newCut

/* This function allocates memory for a new cut structure.  This entails the structure itself, and the _val_ array of oneCut pointers
 * inside the structure.  The actual oneCut structures are allocated according to the numBeta parameter, via calls to new_cut(). */
cutsType *newCuts(int maxCuts) {
	cutsType *cuts;

	if (maxCuts == 0)
		return NULL;

	if (!(cuts = (cutsType *) mem_malloc (sizeof(cutsType))))
		errMsg("allocation", "newCuts", "cuts",0);
	if (!(cuts->vals = (oneCut **) arr_alloc (maxCuts, oneCut)))
		errMsg("allocation", "newCuts", "oneCuts",0);
	cuts->cnt = 0;

	return cuts;
}//END newCuts

/* This function will remove the oldest cut whose corresponding dual variable is zero (thus, a cut which was slack in last solution). */
int reduceCuts(cellType *cell, dVector candidX, dVector pi, int betaLen, double lb) {
	double height, minHeight;
	int minObs, oldestCut,idx;

	minObs 	  = cell->omega->sampleSize;
	oldestCut = cell->cuts->cnt;

	/* identify the oldest loose cut */
	for (idx = 0; idx < cell->cuts->cnt; idx++) {
		if ( idx == cell->iCutIdx || cell->cuts->vals[idx]->rowNum < 0)
			/* avoid dropping incumbent cut and newly added cuts */
			continue;

		if (cell->cuts->vals[idx]->numSamples < minObs && DBL_ABS(pi[cell->cuts->vals[idx]->rowNum + 1]) <= config.TOLERANCE ) {
			minObs = cell->cuts->vals[idx]->numSamples;
			oldestCut = idx;
		}
	}

	/* if the oldest loose cut is the most recently added cut, then the cut with minimium cut height will be dropped */
	if ( oldestCut == cell->cuts->cnt ) {
		minHeight = cutHeight(cell->cuts->vals[0], cell->omega->sampleSize, candidX, betaLen, lb);
		oldestCut = 0;

		for (idx = 1; idx < cell->cuts->cnt; idx++) {
			if (idx == cell->iCutIdx)
				continue;

			height = cutHeight(cell->cuts->vals[idx], cell->omega->sampleSize, candidX, betaLen, lb);
			if (height < minHeight) {
				minHeight = height;
				oldestCut = idx;
			}
		}
	}

	/* drop the selected cut and swap the last cut into its place */
	if ( dropCut(cell, oldestCut) ){
		errMsg("algorithm", "reduceCuts", "failed to drop a cut", 0);
		return -1;
	}

	return oldestCut;
}//END reduceCuts()

/* This function removes a cut from both the cutType structure and the master problem constraint matrix.  In the cuts->vals array, the last
 * cut is swapped into the place of the exiting cut.  In the constraint matrix, the row is deleted, and the row numbers of all constraints
 * below it are decremented. */
int dropCut(cellType *cell, int cutIdx) {
	int idx, deletedRow;

	deletedRow = cell->cuts->vals[cutIdx]->rowNum;
	/* Get rid of the indexed cut on the solver */
	if (  removeRow(cell->master->lp, deletedRow, deletedRow) ) {
		printf("stopped at %d",cell->k);
		errMsg("solver", "dropCut", "failed to remove a row from master problem", 0);
		return 1;
	}
	freeOneCut(cell->cuts->vals[cutIdx]);

	/* move the last cut to the deleted cut's position (structure) */
	cell->cuts->vals[cutIdx] = cell->cuts->vals[--cell->cuts->cnt];

	/* if the swapped cut happens to be the incumbent cut, then update its index */
	if ( cell->iCutIdx == cell->cuts->cnt )
		cell->iCutIdx = cutIdx;

	/* Decrement the row number of all optimality cuts which were added after the cut that was just dropped */
	for (idx = 0; idx < cell->cuts->cnt; idx++) {
		if (cell->cuts->vals[idx]->rowNum > deletedRow)
			--cell->cuts->vals[idx]->rowNum;
	}

	/* Decrement the row number of all feasibility cuts which were added after the cut that was just dropped */
	for (idx = 0; idx < cell->fcuts->cnt; idx++) {
		if (cell->fcuts->vals[idx]->rowNum > deletedRow)
			--cell->fcuts->vals[idx]->rowNum;
	}

	/* decrease the number of rows on solver */
	cell->master->mar--;

	return 0;
}//END dropCut()

/*
 ** This function calculate the variance of the
 ** dVector x.
 */
double calcVariance(double *x, double *mean_value, double *stdev_value, int batch_size) {
	double mean, vari, temp;
	int count, length;
	double stdev;
	stdev = 10000000.0;
	temp = 0.0;
	mean = x[0];
	vari = 0.0;

	if (mean_value != NULL)
		length = batch_size;
	else
		length = config.SCAN_LEN;

	for (count = 1; count < length; count++) {
		temp = mean;
		mean = mean + (x[count] - mean) / (double) (count + 1);
		vari = (1 - 1 / (double) count) * vari
				+ (count + 1) * (mean - temp) * (mean - temp);
	}

	if (mean_value != NULL)
		*mean_value = mean;
	if (stdev_value != NULL) {
		stdev = sqrt(vari / (double) count);
		*stdev_value = stdev;
	}

	return vari;

}//END calcVariance()

/* This function takes the SD code into "feasibility mode, from the "optimality mode". SD will not return to optimality mode
 * until the a feasible candidate and incumbent solution are identified.*/
/* This function takes the SD code into Feasibility mode (solve_cell() take the SD into Optimality mode). The SD will not return to optimality mode
 until the candidate and incumbent solution are both feasible. */
int resolveInfeasibility(probType **prob, cellType *cell, bool *newOmegaFlag, int omegaIdx) {
	bool newBasisFlag;

	/* QP master will be solved in feasibility mode */
	cell->optMode = false;

	while ( true ) {
		/* form a feasibility cut */
		formFeasCut(prob[1], cell);

		/* relax the proximal term and change it in the solver */
		cell->quadScalar = config.MIN_QUAD_SCALAR;
		if ( changeQPproximal(cell->master->lp, prob[0]->num->cols, cell->quadScalar) ) {
			errMsg("algorithm", "resolveInfeasibility", "failed to change the proximal parameter", 0);
			return 1;
		}

		/* Solver the master problem with the added feasibility cut */
		if ( solveQPMaster(prob[0]->num, prob[0]->dBar, cell, prob[0]->lb) ) {
			errMsg("algorithm", "resolveInfeasibility", "failed to solve the master problem", 0);
			return 1;
		}

		/* increment the count for number of infeasible master solutions encountered */
		cell->feasCnt++;

		if ( solveSubprob(prob[1], cell->subprob->lp, cell->candidX, cell->basis, cell->lambda, cell->sigma, cell->delta, config.MAX_ITER,
				cell->omega, omegaIdx, newOmegaFlag, cell->k, config.TOLERANCE, &cell->spFeasFlag, &newBasisFlag,
				&cell->time.subprobIter, &cell->time.argmaxIter) < 0 ) {
			errMsg("algorithm", "resolveInfeasibility", "failed to solve the subproblem", 0);
			return 1;
		}

		/* end the feasibility mode if a feasible candidate solution is observed */
		if (cell->spFeasFlag == true)
			break;
	}

	if ( cell->infeasIncumb == true ) {
		/* if the incumbent solution is infeasible then replace the incumbent with the feasible candidate solution */
		replaceIncumbent(prob[0], cell, cell->candidEst);
	}

	/* QP master will be solved in optimality mode again */
	cell->optMode = true;

	return 0;
}//END resolveInfeasibility()

int formFeasCut(probType *prob, cellType *cell) {

	/* add new feasibility cuts to the cut pool */
	updtFeasCutPool(prob->num, prob->coord, cell);

	/* identify, in the feasibility cut pool, cuts that are violated by the input solution xk */
	checkFeasCutPool(cell, prob->num->prevCols);
	cell->piM = (dVector) mem_realloc(cell->piM, prob->num->prevRows+cell->cuts->cnt+cell->fcuts->cnt);

	return 0;
}//END formFeasCut()

/* This function adds new feasibility cuts. It first adds feasibility cuts from old pi's associated with the new omega generated.
 * Cuts from a new dual extreme ray(new pi) and all omegas generated so far are added to the feasible_cuts_pool structure afterwards. */
int updtFeasCutPool(numType *num, coordType *coord, cellType *cell) {
	oneCut	*cut;
	int		idx, obs, c, initCutsCnt, sigmaIdx, lambdaIdx;

	initCutsCnt = cell->fcutsPool->cnt;

	/* Update computations with respect to the newly discovered observations and all the elements of the stochastic structures. */
	for ( obs = cell->fUpdt[1]; obs < cell->omega->cnt; obs++ )
		for ( idx = 0; idx < cell->fUpdt[0]; idx++ ) {
			if ( !cell->basis->vals[idx]->feasFlag ) {
				cut = newCut(num->prevCols, 0, 1);

				sigmaIdx = cell->basis->vals[idx]->sigmaIdx[0];
				lambdaIdx = cell->sigma->lambdaIdx[sigmaIdx];

				/* Average using these Pi's to calculate the cut itself (update alpha and beta) */
				cut->alpha = cell->sigma->vals[sigmaIdx].pib + cell->delta->vals[lambdaIdx][obs].pib;

				for (c = 1; c <= num->cntCcols; c++)
					cut->beta[coord->CCols[c]] += cell->sigma->vals[sigmaIdx].piC[c];
				for (c = 1; c <= num->rvCOmCnt; c++)
					cut->beta[coord->rvCols[c]] += cell->delta->vals[lambdaIdx][obs].piC[c];

				addCut2Pool(cell, cut, num->prevCols, 0.0, true);
			}
		}
	cell->fUpdt[1] = cell->omega->cnt;

	/* TODO: Update computations with respect to the newly discovered stochastic structures and all the observations discovered until
	 * now. */
	for ( obs = 0; obs < cell->omega->cnt; obs++ )
		for ( idx = cell->fUpdt[0]; idx < cell->basis->cnt; idx++ ) {
			if ( !cell->basis->vals[idx]->feasFlag ) {
				cut = newCut(num->prevCols, 0, 1);

				sigmaIdx = cell->basis->vals[idx]->sigmaIdx[0];
				lambdaIdx = cell->sigma->lambdaIdx[sigmaIdx];

				/* Average using these Pi's to calculate the cut itself (update alpha and beta) */
				cut->alpha = cell->sigma->vals[sigmaIdx].pib + cell->delta->vals[lambdaIdx][obs].pib;

				for (c = 1; c <= num->cntCcols; c++)
					cut->beta[coord->CCols[c]] += cell->sigma->vals[sigmaIdx].piC[c];
				for (c = 1; c <= num->rvCOmCnt; c++)
					cut->beta[coord->rvCols[c]] += cell->delta->vals[lambdaIdx][obs].piC[c];

				addCut2Pool(cell, cut, num->prevCols, 0.0, true);
			}
		}
	cell->fUpdt[0] = cell->basis->cnt;

	return (cell->fcutsPool->cnt - initCutsCnt);
}//END updtFeasCutPool()

/* The function identifies cuts from the feasibility cut pool which are voilated by the candidate solution, and mark them to be
 * added to master problem. */
int checkFeasCutPool(cellType *cell, int lenX) {
	double 	betaX, alpha;
	int 	idx, c;
	bool 	duplicCut;

	for (idx = 0; idx < cell->fcutsPool->cnt; idx++) {
		duplicCut = false;
		alpha = cell->fcutsPool->vals[idx]->alpha;
		for (c = 0; c < cell->fcuts->cnt; c++) {
			if (DBL_ABS(alpha - cell->fcuts->vals[c]->alpha) < config.TOLERANCE) {
				if (equalVector(cell->fcutsPool->vals[idx]->beta, cell->fcuts->vals[c]->beta, lenX, config.TOLERANCE)) {
					duplicCut = true;
					break;
				}
			}
		}

		/* Add those cuts in cut pool that will be violated by incumbent solution */
		betaX = vXv(cell->fcutsPool->vals[idx]->beta, cell->incumbX, NULL, lenX);
		if (betaX < alpha) {
			cell->infeasIncumb = true;
			if (duplicCut == true)
				printf("Incumbent violates one old cut from feasible cut pool, but the cut already exists in the master problem.\n");
			else {
				printf( "Incumbent violates one new cut from feasible cut pool, also adding to master problem.\n");
				addCut2Master(cell->master, cell->fcutsPool->vals[idx], cell->incumbX, lenX);
			}
		}
		else {
			/* Check if the cut will be violated by the candidate solution*/
			if (duplicCut == true)
				continue;
			betaX = vXv(cell->fcutsPool->vals[idx]->beta, cell->candidX, NULL, lenX);

			if (betaX < alpha) {
				printf("Candidate violates one cut from feasible cut pool (this cut is not in feasCutsAdded but will be added)\n");
				addCut2Master(cell->master, cell->fcutsPool->vals[idx], cell->incumbX, lenX);
			}
		}
	}

#if defined(ALGO_CHECK)
	writeProblem(cell->master->lp, "cellMasterAfterFeasCuts.lp");
#endif

	return 0;
}//END checkFeasCutPool()

/* This function prints the relevant information in a cut. It is meant to be used for debugging. */
void printCut(oneCut *cut, int betaLen) {

	printf("Sample size = %d; Obs = %d; Alpha = %.3lf; Beta = ", cut->numSamples, cut->omegaCnt, cut->alpha);
	printVector(cut->beta, betaLen, stdout);
	printf("Istar: ");
	printIntvec(cut->iStar, cut->omegaCnt, stdout);

}//END printCut()


void freeOneCut(oneCut *cut) {

	if (cut) {
		if (cut->iStar)
			mem_free(cut->iStar);
		if (cut->beta)
			mem_free(cut->beta);
		if (cut->name)
			mem_free(cut->name);
		mem_free(cut);
	}
}//END freeOneCut()

void freeCutsType(cutsType *cuts, bool partial) {
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++)
		freeOneCut(cuts->vals[cnt]);

	if ( partial )
		cuts->cnt = 0;
	else {
		mem_free(cuts->vals);
		mem_free(cuts);
	}
}//END freeCutsType()

/* The subroutine adds the newly formed cut to the cutsType structure. For the optimality cuts, if there is no room in the cutsType structure
 * then the reduceCuts() subroutine is invoked to remove the 'loose' and 'old' cuts. For the feasibility cuts, we check if the new cut is a
 * duplicate of existing cut before it is added to the pool. */
int addCut2Pool(cellType *cell, oneCut *cut, int lenX, double lb, bool feasCut) {
	int	cnt;

	if ( feasCut ) {
		/* If we are adding a feasibility cut, make sure there are no duplicates */
		for (cnt = 0; cnt < cell->fcutsPool->cnt; cnt++) {
			if (DBL_ABS(cut->alpha - cell->fcutsPool->vals[cnt]->alpha) < config.TOLERANCE) {
				if (equalVector(cut->beta, cell->fcutsPool->vals[cnt]->beta, lenX, config.TOLERANCE)) {
					/* return 0 to indicate that no cut was added to the pool */
					freeOneCut(cut);
					return 0;
				}
			}
		}
		cell->fcutsPool->vals[cell->fcutsPool->cnt] = cut;
		return cell->fcutsPool->cnt++;
	}
	else {
		if (cell->cuts->cnt >= cell->maxCuts) {
			/* If we are adding optimality cuts, check to see if there is room for the latest cut. If there is not,
			 * then make room by reducing the cuts from the structure. */
			if( reduceCuts(cell, cell->candidX, cell->piM, lenX, lb) < 0 ) {
				errMsg("algorithm", "addCut2Master", "failed to add reduce cuts to make room for candidate cut", 0);
				return -1;
			}
		}
		cell->cuts->vals[cell->cuts->cnt] = cut;
		return cell->cuts->cnt++;
	}

}//END addCut()
