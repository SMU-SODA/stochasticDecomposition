/*
 * cuts.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 *
 *   Edited on: Aug 1, 2020 as part of the SD-integer project
 *      Author: Siavash Tabrizian
 * Institution: Southern Methodist University
 *
 * Please send you comments or bug report to stabrizian (at) smu (dot) edu
 *
 */

#include "twoSD.h"

#undef Disp_argmax  

#undef remove_cut

extern configType config;

int resolveInfeasibility(probType **prob, cellType *cell, bool *newOmegaFlag, int omegaIdx);
int LPresolveInfeasibility(probType **prob, cellType *cell, bool *newOmegaFlag, int omegaIdx);
int formFeasCut(probType *prob, cellType *cell);
int updtFeasCutPool(numType *num, coordType *coord, cellType *cell);
int checkFeasCutPool(cellType *cell, int lenX);
int addCut2ActivePool(cellType *cell, int mar, oneCut *cut, int lenX, double lb, bool feasCut);

int formSDCut(probType **prob, cellType *cell, dVector Xvect, double lb, int inCallback) {
	oneCut 	*cut;
	int    	cutIdx, obs;

	/* A subproblem is solved for every new observation */
	for ( obs = 0; obs < config.SAMPLE_INCREMENT; obs++ ) {
		/* (a) Construct the subproblem with input observation and master solution, solve the subproblem, and complete stochastic updates */
		if ( solveSubprob(prob[1], cell->subprob, Xvect, cell->basis, cell->lambda, cell->sigma, cell->delta,
				config.MAX_ITER, cell->omega, cell->sample->omegaIdx[obs], &cell->sample->newOmegaFlag[obs], cell->k, config.TOLERANCE,
				&cell->spFeasFlag, &cell->sample->newBasisFlag[obs], &cell->time.subprobIter, &cell->time.argmaxIter) < 0 ){
			errMsg("algorithm", "formSDCut", "failed to solve the subproblem", 0);
			return -1;
		}
		else {
			cell->sample->basisIdx[obs] = solveSubprob(prob[1], cell->subprob, Xvect, cell->basis, cell->lambda, cell->sigma, cell->delta,
				config.MAX_ITER, cell->omega, cell->sample->omegaIdx[obs], &cell->sample->newOmegaFlag[obs], cell->k, config.TOLERANCE,
				&cell->spFeasFlag, &cell->sample->newBasisFlag[obs], &cell->time.subprobIter, &cell->time.argmaxIter);
		}

		/* increment the number of subproblems solved during algorithm */
		cell->LPcnt++;

		if ( ! cell->spFeasFlag ) {
			if (inCallback == 1)
			{
				/* Subproblem is infeasible, resolve infeasibility */
				if (LPresolveInfeasibility(prob, cell, &cell->sample->newOmegaFlag[obs], cell->sample->omegaIdx[obs])) {
					errMsg("algorithm", "formSDCut", "failed to resolve infeasibility", 0);
					return -1;
				}
			}
			else
			{
				/* Subproblem is infeasible, resolve infeasibility */
				if (resolveInfeasibility(prob, cell, &cell->sample->newOmegaFlag[obs], cell->sample->omegaIdx[obs])) {
					errMsg("algorithm", "formSDCut", "failed to resolve infeasibility", 0);
					return -1;
				}
			}
		}
		else if ( cell->fcutsPool->cnt > 0 && cell->sample->newOmegaFlag[obs] ) {
			/* Subproblem is feasible, however a new observation or sigma has been encountered. Therefore, update the feasibility cut pool and check
			 * to see if new feasibility cuts need to be added. */
			if ( formFeasCut(prob[1], cell) ) {
				errMsg("algorithm", "formSDCut", "failed to add new feasibility cuts", 0);
				return -1;
			}
		}
	}

	/* (b) create an affine lower bound */
	clock_t tic = clock();
	cut = SDCut(prob[1]->num, prob[1]->coord, cell->basis, cell->sigma, cell->delta, cell->omega, cell->sample,
			Xvect, cell->sampleSize, &cell->dualStableFlag, cell->pi_ratio, cell->k, cell->lb);
	if ( cut == NULL ) {
		errMsg("algorithm", "formSDCut", "failed to create the affine minorant", 0);
		return -1;
	}
	cell->time.argmaxIter += ((double) (clock()-tic))/CLOCKS_PER_SEC;

#if defined(STOCH_CHECK)
	/* Solve the subproblem to verify if the argmax operation yields a lower bound */
	for ( int cnt = 0; cnt < cell->omega->cnt; cnt++ ) {
		/* (a) Construct the subproblem with input observation and master solution, solve the subproblem, and complete stochastic updates */
		if ( solveSubprob(prob[1], cell->subprob, Xvect, cell->basis, cell->lambda, cell->sigma, cell->delta, config.MAX_ITER,
				cell->omega, cnt, newOmegaFlag, cell->k, config.TOLERANCE, &cell->spFeasFlag, NULL,
				&cell->time.subprobIter, &cell->time.argmaxIter) < 0 ) {
			errMsg("algorithm", "formSDCut", "failed to solve the subproblem", 0);
			return -1;
		}
		printf("Subproblem solve for omega-%d = %lf\n", cnt, getObjective(cell->subprob->lp, PROB_LP));
	}
#endif

	/* (c) add cut to the structure and master problem  */
	if ( (cutIdx = addCut2ActivePool(cell, prob[0]->num->rows, cut, prob[0]->num->cols, lb, false)) < 0) {
		errMsg("algorithm", "formSDCut", "failed to add the new cut to cutsType structure", 0);
		return -1;
	}

	if ( addCut2Master(cell->master, cut, cell->incumbX, prob[0]->num->cols, cell->master->type == PROB_QP) ) {
		errMsg("algorithm", "formSDCut", "failed to add the new cut to master problem", 0);
		return -1;
	}

	return cutIdx;
}//END formCut()

oneCut *SDCut(numType *num, coordType *coord, basisType *basis, sigmaType *sigma, deltaType *delta, omegaType *omega, sampleType *sample,
		dVector Xvect, int numSamples, bool *dualStableFlag, dVector pi_ratio, int numIter, double lb) {
	oneCut *cut;
	dVector 	piCbarX, beta;
	double  argmaxOld, argmaxNew, cummOld = 0.0, cummAll = 0.0, argmax, alpha = 0.0, variance = 1.0, multiplier;
	int	 	istarOld, istarNew, istar, idx, c, obs, sigmaIdx, lambdaIdx;
	bool    pi_eval_flag = false;

	/* allocate memory to hold a new cut */
	cut = newCut(num->prevCols, omega->cnt, numSamples);

	/* Pre-compute pi x Cbar x x as it is independent of observations */
	if (!(piCbarX= arr_alloc(sigma->cnt, double)))
		errMsg("Allocation", "SDCut", "pi_Tbar_x",0);
	for (c = 0; c < sigma->cnt; c++)
		piCbarX[c] = vXv(sigma->vals[c].piC, Xvect, coord->CCols, num->cntCcols);

	if ( !(beta = (dVector) arr_alloc(num->prevCols + 1, double)) )
		errMsg("Allocation", "SDCut", "beta", 0);

	/* Calculate pi_eval_flag to determine the way of computing argmax */
	if (config.DUAL_STABILITY && numSamples > config.PI_EVAL_START && !(numSamples % config.PI_CYCLE))
		pi_eval_flag = true;

	/* Test for omega issues */
	for (obs = 0; obs < omega->cnt; obs++) {
		/* For each observation, find the Pi/basis that generates the Pi which maximizes height at X. */
		if (pi_eval_flag == true) {
			istarOld = computeIstar(num, coord, basis, sigma, delta, sample,
					piCbarX, Xvect, omega->vals[obs], obs, numSamples, pi_eval_flag, &argmaxOld, false);
			istarNew = computeIstar(num, coord, basis, sigma, delta, sample,
					piCbarX, Xvect, omega->vals[obs], obs, numSamples, true, &argmaxNew, true);

			argmax = maximum(argmaxOld, argmaxNew);

#if defined(Disp_argmax)
	printf("\nargmax val: %0.4f", argmax);
#endif // define(Disp_argmax)
			
			istar  = (argmaxNew > argmaxOld) ? istarNew : istarOld;

			cummOld += maximum(argmaxOld-lb, 0)*omega->weights[obs];
			cummAll += maximum(argmax-lb, 0)*omega->weights[obs];
		}
		else {
			/* identify the maximal Pi/basis that generates the maximal Pi for each observation */
			istar = computeIstar(num, coord, basis, sigma, delta, sample,
					piCbarX, Xvect, omega->vals[obs], obs, numSamples, pi_eval_flag, &argmax, false);
		}

		if (istar < 0) {
			errMsg("algorithm", "SDCut", "failed to identify maximal Pi for an observation", 0);
			return NULL;
		}
		cut->iStar[obs] = istar;

		if ( num->rvdOmCnt > 0 ) {
			for ( idx = 0; idx <= basis->vals[istar]->phiLength; idx++ ) {
				sigmaIdx = basis->vals[istar]->sigmaIdx[idx];
				lambdaIdx = sigma->lambdaIdx[sigmaIdx];
				if ( idx == 0 )
					multiplier = 1.0;
				else
					multiplier = omega->vals[obs][coord->rvOffset[2] + basis->vals[istar]->omegaIdx[idx]];

				/* Start with (Pi x bBar) + (Pi x bomega) + (Pi x Cbar) x X */
				alpha += omega->weights[obs] * multiplier * (sigma->vals[sigmaIdx].pib + delta->vals[lambdaIdx][obs].pib);

				for (c = 1; c <= num->cntCcols; c++)
					beta[coord->CCols[c]] += omega->weights[obs] * multiplier * sigma->vals[sigmaIdx].piC[c];
				for (c = 1; c <= num->rvCOmCnt; c++)
					beta[coord->rvCOmCols[c]] += omega->weights[obs] * multiplier * delta->vals[lambdaIdx][obs].piC[c];
			}
		}
		else {
			alpha += sigma->vals[istar].pib * omega->weights[obs];
			alpha += delta->vals[sigma->lambdaIdx[istar]][obs].pib * omega->weights[obs];

			for (c = 1; c <= num->cntCcols; c++)
				beta[coord->CCols[c]] += sigma->vals[istar].piC[c] * omega->weights[obs];
			for (c = 1; c <= num->rvCOmCnt; c++)
				beta[coord->rvCols[c]] += delta->vals[sigma->lambdaIdx[istar]][obs].piC[c] * omega->weights[obs];
		}
	}

	if (pi_eval_flag == true) {
		pi_ratio[numIter % config.SCAN_LEN] = cummOld / cummAll;
		if (numSamples - config.PI_EVAL_START > config.SCAN_LEN)
			variance = calcVariance(pi_ratio, NULL, NULL, 0);
		else
			variance = 1.0;

		if (DBL_ABS(variance) >= .000002 || pi_ratio[numSamples % config.SCAN_LEN] < 0.95 )
			*dualStableFlag = false;
		else
			*dualStableFlag = true;
	}

	cut->alpha = alpha / numSamples;

	for (c = 1; c <= num->prevCols; c++)
		cut->beta[c] = beta[c] / numSamples;
	cut->beta[0] = 1.0;			/* coefficient of eta coloumn - it is the last column in the master problem but the first element 
								in beta vector, should look at add2Master subroutine */

	mem_free(piCbarX);
	mem_free(beta);

	return cut;
}//END SDCut


/* This function loops through a set of cuts and find the highest cut height at the specified position x */
double maxCutHeight(cutsType *cuts, int currSampleSize, dVector xk, int betaLen, double lb) {
	double Sm = -INF, ht = 0.0;
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		if (cuts->vals[cnt])
		{
			ht = cutHeight(cuts->vals[cnt], currSampleSize, xk, betaLen, lb);
			if (Sm < ht) {
				Sm = ht;
			}
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

	if (!(cut->iStar = arr_alloc(config.MAX_ITER, int)))
		errMsg("allocation", "new_cut", "iStar", 0);

	if (!(cut->beta = arr_alloc(numX + 1, double)))
		errMsg("allocation", "new_cut", "beta", 0);

	cut->alpha = 0.0;

	cut->name = (cString) arr_alloc(3*NAMESIZE, char);

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

	minObs 	  = cell->sampleSize;
	oldestCut = cell->activeCuts->cnt;

	/* identify the oldest loose cut */
	for (idx = 0; idx < cell->activeCuts->cnt; idx++) {
		if ( idx == cell->iCutIdx || cell->activeCuts->vals[idx]->rowNum < 0)
			/* avoid dropping incumbent cut and newly added cuts */
			continue;

		if (cell->activeCuts->vals[idx]->numSamples < minObs && DBL_ABS(pi[cell->activeCuts->vals[idx]->rowNum + 1]) <= config.TOLERANCE ) {
			minObs = cell->activeCuts->vals[idx]->numSamples;
			oldestCut = idx;
		}
	}

	/* if the oldest loose cut is the most recently added cut, then the cut with minimium cut height will be dropped */
	if ( oldestCut == cell->activeCuts->cnt ) {
		minHeight = cutHeight(cell->activeCuts->vals[0], cell->sampleSize, candidX, betaLen, lb);
		oldestCut = 0;

		for (idx = 1; idx < cell->activeCuts->cnt; idx++) {
			if (idx == cell->iCutIdx)
				continue;

			height = cutHeight(cell->activeCuts->vals[idx], cell->sampleSize, candidX, betaLen, lb);
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

	deletedRow = cell->activeCuts->vals[cutIdx]->rowNum;
	/* Get rid of the indexed cut on the solver */
	if (  removeRow(cell->master->lp, deletedRow, deletedRow) ) {
		printf("stopped at %d",cell->k);
		errMsg("solver", "dropCut", "failed to remove a row from master problem", 0);
		return 1;
	}
	freeOneCut(cell->activeCuts->vals[cutIdx]);

	/* move the last cut to the deleted cut's position (structure) */
	cell->activeCuts->vals[cutIdx] = cell->activeCuts->vals[--cell->activeCuts->cnt];

	/* if the swapped cut happens to be the incumbent cut, then update its index */
	if ( cell->iCutIdx == cell->activeCuts->cnt )
		cell->iCutIdx = cutIdx;

	/* Decrement the row number of all optimality cuts which were added after the cut that was just dropped */
	for (idx = 0; idx < cell->activeCuts->cnt; idx++) {
		if (cell->activeCuts->vals[idx]->rowNum > deletedRow)
			--cell->activeCuts->vals[idx]->rowNum;
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

int copyCuts(numType *num, cutsType *orig, cutsType **copy) {


	if ( (*copy) == NULL ) {
		(*copy) = newCuts(orig->cnt);
	}

	for ( int cnt = 0; cnt < orig->cnt; cnt++ ) {
		oneCut *cut;

		cut = newCut(num->cols, orig->vals[cnt]->omegaCnt, orig->vals[cnt]->numSamples);
		copyOneCut(orig->vals[cnt], cut, num->cols);

		(*copy)->vals[(*copy)->cnt++] = cut;
	}

	(*copy)->cnt = orig->cnt;

	return 0;
}//END copyCuts()

void copyOneCut(oneCut *orig, oneCut *copy, int numCols) {

	copy->isIncumb = orig->isIncumb;
	copy->alpha = orig->alpha;
	copy->alphaIncumb = orig->alphaIncumb;
	copy->slackCnt = orig->slackCnt;

	copyVector(orig->beta, copy->beta, numCols, true);
	copyIntvec(orig->iStar, copy->iStar, orig->omegaCnt);
	strcpy(copy->name, orig->name);

}//END copyOneCut()

int copyCutstoNodePool(numType *num, cutsType *orig, cutsType *copy, dVector pi) {

	/* If the current pool is filled, clean it before refilling. This happens as one of the daughter node inherits the
	 * parent's pool. */
	if ( copy->cnt > 0 ) {
		freeCutsType(copy, true);
	}

	int count = 0;
	if (orig != NULL && orig->cnt > 0) { /* if the parent cuts are not empty continue */
		for (int cnt = 0; cnt < orig->cnt; cnt++) {
			if (pi[orig->vals[cnt]->rowNum + 1] > config.TOLERANCE) {
				oneCut *cut;

				cut = (oneCut *)mem_malloc(sizeof(oneCut));
				cut->numSamples = orig->vals[cnt]->numSamples;
				cut->omegaCnt = orig->vals[cnt]->omegaCnt;

				cut->rowNum = -1;
				cut->isIncumb = orig->vals[cnt]->isIncumb;
				cut->alpha = orig->vals[cnt]->alpha;
				cut->alphaIncumb = orig->vals[cnt]->alphaIncumb;
				cut->slackCnt = orig->vals[cnt]->slackCnt;

				cut->beta = duplicVector(orig->vals[cnt]->beta, num->cols + 1, false);
				if (orig->vals[cnt]->iStar == NULL) {
					cut->iStar = NULL;
				}
				else {
					cut->iStar = duplicIntvec(orig->vals[cnt]->iStar, orig->vals[cnt]->omegaCnt, false);
				}

				cut->name = (cString)arr_alloc(NAMESIZE, char);
				strcpy(cut->name, orig->vals[cnt]->name);

				copy->vals[copy->cnt++] = cut;

				count++;
			}
		}

		if (count == 0) {
			oneCut *cut;

			cut = (oneCut *)mem_malloc(sizeof(oneCut));
			cut->numSamples = orig->vals[0]->numSamples;
			cut->omegaCnt = orig->vals[0]->omegaCnt;

			cut->rowNum = -1;
			cut->isIncumb = orig->vals[0]->isIncumb;
			cut->alpha = orig->vals[0]->alpha;
			cut->alphaIncumb = orig->vals[0]->alphaIncumb;
			cut->slackCnt = orig->vals[0]->slackCnt;

			cut->beta = duplicVector(orig->vals[0]->beta, num->cols + 1, false);
			if (orig->vals[0]->iStar == NULL) {
				cut->iStar = NULL;
			}
			else {
				cut->iStar = duplicIntvec(orig->vals[0]->iStar, orig->vals[0]->omegaCnt, false);
			}

			cut->name = (cString)arr_alloc(NAMESIZE, char);
			strcpy(cut->name, orig->vals[0]->name);

			copy->vals[copy->cnt++] = cut;

		}
	}
	else {
		errMsg("cleanNode", "copyCutstoNodePool", "orig cut pool is empty", 0);
	}

	return 0;
}//END copyCutstoNodePool()

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

/* resolve infeasibility for LP relaxation 
   added by Siavash Tabrizian July 20
*/
int LPresolveInfeasibility(probType **prob, cellType *cell, bool *newOmegaFlag, int omegaIdx) {
	bool newBasisFlag;

	/* LP master will be solved in feasibility mode */
	cell->optMode = false;

	while (true) {
		/* form a feasibility cut */
		formFeasCut(prob[1], cell);

		
		/* Solver the master problem with the added feasibility cut */
		if (solveLPMaster(prob[0]->num, prob[0]->dBar, cell, prob[0]->lb)) {
			errMsg("algorithm", "resolveInfeasibility", "failed to solve the master problem", 0);
			return 1;
		}

		/* increment the count for number of infeasible master solutions encountered */
		cell->feasCnt++;

		if (solveSubprob(prob[1], cell->subprob->lp, cell->candidX, cell->basis, cell->lambda, cell->sigma, cell->delta, config.MAX_ITER,
			cell->omega, omegaIdx, newOmegaFlag, cell->k, config.TOLERANCE, &cell->spFeasFlag, &newBasisFlag,
			&cell->time.subprobIter, &cell->time.argmaxIter) < 0) {
			errMsg("algorithm", "resolveInfeasibility", "failed to solve the subproblem", 0);
			return 1;
		}

		/* end the feasibility mode if a feasible candidate solution is observed */
		if (cell->spFeasFlag == true)
			break;
	}

	if (cell->infeasIncumb == true) {
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
	cell->piM = (dVector) mem_realloc(cell->piM, prob->num->prevRows + cell->activeCuts->cnt+cell->fcuts->cnt);

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

				addCut2ActivePool(cell, 0, cut, num->prevCols, 0.0, true);
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

				addCut2ActivePool(cell, 0, cut, num->prevCols, 0.0, true);
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
				addCut2Master(cell->master, cell->fcutsPool->vals[idx], cell->incumbX, lenX, cell->master->type == PROB_QP);
			}
		}
		else {
			/* Check if the cut will be violated by the candidate solution*/
			if (duplicCut == true)
				continue;
			betaX = vXv(cell->fcutsPool->vals[idx]->beta, cell->candidX, NULL, lenX);

			if (betaX < alpha) {
				printf("Candidate violates one cut from feasible cut pool (this cut is not in feasCutsAdded but will be added)\n");
				addCut2Master(cell->master, cell->fcutsPool->vals[idx], cell->incumbX, lenX, cell->master->type == PROB_QP);
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
int addCut2ActivePool(cellType *cell, int mar, oneCut *cut, int lenX, double lb, bool feasCut) {
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
		if (cell->master->mar - mar >= cell->maxCuts) {
			/* If we are adding optimality cuts, check to see if there is room for the latest cut. If there is not,
			 * then make room by reducing the cuts from the structure. */
			if( reduceCuts(cell, cell->candidX, cell->piM, lenX, lb) < 0 ) {
				errMsg("algorithm", "addCut2Master", "failed to add reduce cuts to make room for candidate cut", 0);
				return -1;
			}
		}

		cell->activeCuts->vals[cell->activeCuts->cnt] = cut;
		return cell->activeCuts->cnt++;
	}

}//END addCut()
