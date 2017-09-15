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

#include "cell.h"

extern configType config;

int addCut2Master(cellType *cell, oneCut *cut, BOOL scaleCut, int lenX, double lb) {
	intvec 	indices;
	int 	cnt;

	if (!(indices = arr_alloc(lenX + 1, int)))
		errMsg("Allocation", "addcut2Master", "fail to allocate memory to coefficients of beta",0);
	for (cnt = 1; cnt <= lenX; cnt++)
		indices[cnt] = cnt - 1;
	indices[0] = lenX;

	if ( config.MASTER_TYPE == PROB_QP )
		cut->alphaIncumb = cut->alpha - vXv(cut->beta, cell->incumbX, NULL, lenX);

	/* This is an optimality cut being added */
	/* check to see if there is room for the candidate cut, else drop a cut */
	if (cell->cuts->cnt == cell->maxCuts) {
		/* make room for the latest cut */
		if( reduceCuts(cell->master, cell->cuts, scaleCut, cell->candidX, cell->piM, lenX, lb, cell->k, &cell->iCutIdx, config.TOLERANCE) < 0 ) {
			errMsg("algorithm", "addCut2Master", "failed to add reduce cuts to make room for candidate cut", 0);
			return -1;
		}
	}

	/* Add the cut to the cell cuts structure and assign a row number. */
	cell->cuts->vals[cell->cuts->cnt] = cut;
	cut->rowNum = cell->master->mar++;

	/* Add the row in the solver */
	if ( addRow(cell->master->lp, lenX + 1, cut->alphaIncumb, GE, 0, indices, cut->beta) ) {
		errMsg("solver", "addcut2Master", "failed to add new row to problem in solver", 0);
		return -1;
	}

	mem_free(indices);
	return cell->cuts->cnt++;
}//END addCuts2Master()

int replaceIncumbent(probType *prob, cellType *cell, double candidEst) {

	/* replace the incumbent solution with the candidate solution */
	copyVector(cell->candidX, cell->incumbX, prob->num->cols, 1);
	cell->incumbEst = candidEst;

	/* update the proximal parameter based on estimated improvement */
	if ( cell->normDk > config.TOLERANCE )
		if ( cell->normDk >= config.R3 * cell->normDk_1 ) {
			cell->quadScalar *= config.R2 * config.R3 * cell->normDk_1/ cell->normDk;
			cell->quadScalar  = min(config.MAX_QUAD_SCALAR, cell->quadScalar);
			cell->quadScalar = max(config.MIN_QUAD_SCALAR, cell->quadScalar);
		}

	/* update the right-hand side and the bounds with new incumbent solution */
	if ( constructQP(prob, cell, cell->incumbX, cell->quadScalar) ) {
		errMsg("algorithm", "replaceIncumbent", "failed to change the right-hand side after incumbent change", 0);
		return 1;
	}

	/* update the candidate cut as the new incumbent cut */
	cell->iCutUpdt = cell->k;
	cell->incumbChg = TRUE;

	/* keep the two norm of solution*/
	cell->normDk_1 = cell->normDk;
	/* Since incumbent solution is now replaced by a candidate, we assume it is feasible now */
	cell->infeasIncumb = FALSE;
	/* gamma needs to be reset to 0 since there's no difference between candidate and incumbent*/
	cell->gamma = 0.0;

	return 0;
}//END replaceIncumbent()

/* This function loops through a set of cuts and find the highest cut height at the specified position x */
double maxCutHeight(cutsType *cuts, vector xk, int betaLen, BOOL scaleCut, int currIter, double lb) {
	double Sm = -INF, ht = 0.0;
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		ht = cutHeight(cuts->vals[cnt], xk, betaLen, scaleCut, currIter, lb);
		if (Sm < ht) {
			Sm = ht;
		}
	}

	return Sm;
}//END maxCutHeight

/* This function calculates and returns the height of a given cut at a given X.  It includes the k/(k-1) update, but does not include
 * the coefficients due to the cell. */
double cutHeight(oneCut *cut, vector xk, int betaLen, BOOL scaleCut, int currIter, double lb) {
	double height;

	/* A cut is calculated as alpha - beta x X */
	height = cut->alpha - vXv(cut->beta, xk, NULL, betaLen);

	if ( scaleCut ) {
		double t_over_k = ((double) cut->numSamples / (double) currIter);

		/* Weight cut based on number of observations used to form it */
		height *= t_over_k;

		/* Updated for optimality cut height*/
		height += (1 - t_over_k) * lb;
	}
	return height;
}//END cutHeight()

/* This function allocates memory for the arrays inside a single cut, and initializes its values accordingly.  The cut structure
 * itself is assumed to be already allocated.  Note, each beta vector contains room for its one-norm, thought it just gets filled
 * with zero anyway. */
oneCut *newCut(int numX, int numIstar, int numSamples) {
	oneCut *cut;

	cut = (oneCut *) mem_malloc (sizeof(oneCut));
	cut->numSamples = numSamples;
	cut->omegaCnt = numIstar;
	cut->isIncumb = FALSE; 								/* new cut is by default not an incumbent */
	cut->alphaIncumb = 0.0;
	cut->rowNum = -1;

	if (!(cut->iStar = arr_alloc(numIstar, int)))		/* when used in aggregate cut mode (MULTI_CUT = 0), this holds the index of agent cuts */
		errMsg("allocation", "new_cut", "iStar", 0);
	if (!(cut->beta = arr_alloc(numX + 1, double)))
		errMsg("allocation", "new_cut", "beta", 0);

	cut->alpha = 0.0;

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
int reduceCuts(oneProblem *master, cutsType *cuts, BOOL scaleCut, vector candidX, vector pi, int betaLen, double lb, int currentIter, int *iCutIdx, double TOLERANCE) {
	double height, minHeight;
	int minObs, oldestCut,idx;

	minObs 	  = currentIter;
	oldestCut = cuts->cnt;

	/* identify the oldest loose cut */
	for (idx = 0; idx < cuts->cnt; idx++) {
		if ( idx == (*iCutIdx) || cuts->vals[idx]->rowNum < 0)
			/* avoid dropping incumbent cut and newly added cuts */
			continue;

		if (cuts->vals[idx]->numSamples < minObs && DBL_ABS(pi[cuts->vals[idx]->rowNum + 1]) <= TOLERANCE ) {
			minObs = cuts->vals[idx]->numSamples;
			oldestCut = idx;
		}
	}

	/* if the oldest loose cut is the most recently added cut, then the cut with minimium cut height will be dropped */
	if ( oldestCut == cuts->cnt ) {
		//minHeight = cutHeight(lbType, cell[agentIdx]->cuts->vals[0], cell[agentIdx]->k, candidX, betaLen, lb);
		minHeight = cutHeight(cuts->vals[0], candidX, betaLen, scaleCut, currentIter, lb);
		oldestCut = 0;

		for (idx = 1; idx < cuts->cnt; idx++) {
			if (idx == (*iCutIdx))
				continue;

			//height = cutHeight(lbType, cell[agentIdx]->cuts->vals[idx], cell[agentIdx]->k, candidX, betaLen, lb);
			height = cutHeight(cuts->vals[idx], candidX, betaLen, scaleCut, currentIter, lb);
			if (height < minHeight) {
				minHeight = height;
				oldestCut = idx;
			}
		}
	}

	/* drop the selected cut and swap the last cut into its place */
	if ( dropCut(master, cuts, oldestCut, iCutIdx) ){
		errMsg("algorithm", "reduceCuts", "failed to drop a cut", 0);
		return -1;
	}

	return oldestCut;
}//END reduceCuts()

/* This function removes a cut from both the cutType structure and the master problem constraint matrix.  In the cuts->vals array, the last
 * cut is swapped into the place of the exiting cut.  In the constraint matrix, the row is deleted, and the row numbers of all constraints
 * below it are decremented. */
int dropCut(oneProblem *master, cutsType *cuts, int cutIdx, int *iCutIdx) {
	int idx, deletedRow;

	deletedRow = cuts->vals[cutIdx]->rowNum;
	/* Get rid of the indexed cut on the solver */
	if (  removeRow(master->lp, deletedRow, deletedRow) ) {
		errMsg("solver", "dropCut", "failed to remove a row from master problem", 0);
		return 1;
	}
	freeOneCut(cuts->vals[cutIdx]);

	/* move the last cut to the deleted cut's position (structure) */
	cuts->vals[cutIdx] = cuts->vals[--cuts->cnt];

	/* if the swapped cut happens to be the incumbent cut, then update its index */
	if ( (*iCutIdx) == cuts->cnt )
		(*iCutIdx) = cutIdx;

	for (idx = 0; idx < cuts->cnt; idx++) {
		if (cuts->vals[idx]->rowNum > deletedRow)
			--cuts->vals[idx]->rowNum;
	}

	/* decrease the number of rows on solver */
	master->mar--;

	return 0;
}//END dropCut()

int formFeasCut(probType *prob, cellType *cell, BOOL *newOmegaFlag, BOOL newBasisFlag) {
	int start, end;
	int idx;

	/* add new feasibility cuts to the cut pool */
	updtFeasCutPool(prob->num, prob->coord, cell->fCutsPool, cell->fUpdt, cell->basis, cell->sigma, cell->delta, cell->omega,
			(*newOmegaFlag), newBasisFlag, cell->k);

	/* identify, in the feasibility cut pool, cuts that are violated by the input solution xk */
	start = cell->fCuts->cnt;
	checkFeasCutPool(cell->fCutsPool, cell->fCuts, prob->num->prevCols, cell->incumbX, cell->candidX, &cell->infeasIncumb);
	end = cell->fCuts->cnt;

	/* add feasibility cuts to master problem */
	if (end > start) {
		for (idx = start; idx < end; idx++) {
			addfCut2Master(cell->master->lp, cell->fCuts->vals[idx], cell->incumbX, prob->num->prevCols, cell->cuts->cnt, idx);
			writeProblem(cell->master->lp, "feasMaster.lp");
		}
		/* make room for dual solutions for new feasibility cuts added */
		cell->piM = (vector) mem_realloc(cell->piM, prob->num->prevRows+cell->cuts->cnt+cell->fCuts->cnt);
	}

	return 0;
}//END formFeasCut()

/*********************************************************************************************
 This function adds new feasibility cuts. It first adds feasibility cuts from old pi's
 associated with the new omega generated. Cuts from a new dual extreme ray(new pi) and all omegas
 generated so far are added to the feasible_cuts_pool structure afterwards.
 *********************************************************************************************/
int updtFeasCutPool(numType *num, coordType *coord, cutsType *fCutsPool, int fUpdt[2], basisType *basis, sigmaType *sigma, deltaType *delta, omegaType *omega,
		BOOL newOmegaFlag, BOOL newBasisFlag, int currentIter) {
	vector 	beta;
	double	alpha, multiplier;
	int		base, obs, cnt, c, cutCnt = 0, offset;

	offset = num->rvbOmCnt + num->rvCOmCnt;
	if ( newOmegaFlag ) {
		for ( obs = fUpdt[1]; obs < omega->cnt; obs++ )
			for ( base = 0; base < fUpdt[0]; base++ ) {
				if ( !(beta = (vector) arr_alloc(num->prevCols+1, double)) )
					errMsg("allocation", "updtFeasCutPool", "beta", 0);
				for ( cnt = 0; cnt <= basis->vals[base]->phiLength; cnt++ ) {
					if (cnt == 0 )
						multiplier = 1.0;
					else
						multiplier = omega->vals[obs][offset+basis->vals[base]->omegaIdx[cnt]];

					/* Average using these Pi's to calculate the cut itself (update alpha and beta) */
					alpha = (sigma->vals[basis->vals[base]->sigmaIdx[cnt]].pib + delta->vals[basis->vals[base]->lambdaIdx[cnt]][obs].pib)* multiplier;

					for (c = 1; c <= num->cntCcols; c++)
						beta[coord->colsC[c]] += sigma->vals[basis->vals[base]->sigmaIdx[cnt]].piC[c] * multiplier;
					for (c = 1; c <= num->rvCOmCnt; c++)
						beta[coord->rvCols[c]] += delta->vals[basis->vals[base]->lambdaIdx[cnt]][obs].piC[c] * multiplier;
				}
				cutCnt += add2CutPool(fCutsPool, alpha, beta, num->prevCols, omega->cnt, currentIter);
			}
		fUpdt[1] = omega->cnt;
	}

	if ( newBasisFlag ) {
		for ( obs = 0; obs < omega->cnt; obs++ )
			for ( base = fUpdt[0]; base < basis->cnt; base++ ) {
				if ( !(beta = (vector) arr_alloc(num->prevCols+1, double)) )
					errMsg("allocation", "updtFeasCutPool", "beta", 0);
				for ( cnt = 0; cnt <= basis->vals[base]->phiLength; cnt++ ) {
					if (cnt == 0 )
						multiplier = 1.0;
					else
						multiplier = omega->vals[obs][offset+basis->vals[base]->omegaIdx[cnt]];

					/* Average using these Pi's to calculate the cut itself (update alpha and beta) */
					alpha = (sigma->vals[basis->vals[base]->sigmaIdx[cnt]].pib + delta->vals[basis->vals[base]->lambdaIdx[cnt]][obs].pib)* multiplier;

					for (c = 1; c <= num->cntCcols; c++)
						beta[coord->colsC[c]] += sigma->vals[basis->vals[base]->sigmaIdx[cnt]].piC[c] * multiplier;
					for (c = 1; c <= num->rvCOmCnt; c++)
						beta[coord->rvCols[c]] += delta->vals[basis->vals[base]->lambdaIdx[cnt]][obs].piC[c] * multiplier;
				}
				cutCnt += add2CutPool(fCutsPool, alpha, beta, num->prevCols, omega->cnt, currentIter);
			}
		fUpdt[0] = basis->cnt;
	}

	return cutCnt;
}//END updtFeasCutPool()

/* This function add a new cut to the cut pool using alpha and beta provided. */
int add2CutPool(cutsType *cuts, double alpha, vector beta, int betaLen, int numOmega, int numSamples) {
	oneCut 	*cut;
	int 	cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		if (DBL_ABS(alpha - cuts->vals[cnt]->alpha) < config.TOLERANCE) {
			if (equalVector(beta, cuts->vals[cnt]->beta, betaLen, config.TOLERANCE)) {
				/* return 0 to indicate that no cut was added to the pool */
				mem_free(beta);
				return 0;
			}
		}
	}

	if ( !(cut = (oneCut *) mem_malloc (sizeof(oneCut))))
		errMsg("allocation", "add2CutPool", "cut", 0);
	cut->numSamples = numSamples;
	cut->omegaCnt = numOmega;
	cut->isIncumb = FALSE;

	if ( !(cut->iStar = (intvec) arr_alloc(numOmega, int)) )
		errMsg("allocation", "add2CutPool", "istar", 0);
	if ( !(cut->beta = arr_alloc(betaLen+1, double)))
		errMsg("allocation", "add2CutPool", "beta", 0);

	cut->alpha = alpha;
	cut->beta = beta;

	cuts->vals[cuts->cnt++] = cut;

	return 1;
}//END add2CutPool()


/* The function identifies cuts from the feasibility cut pool which are voilated by the candidate solution, and mark them to be
 * added to master problem. */
int checkFeasCutPool(cutsType *cutPool, cutsType *cutsAdded, int betaLen, vector incumbX, vector candidX, BOOL *infeasIncumb) {
	double 	betaX, alpha;
	int 	idx, c;
	BOOL 	duplicCut;

	for (idx = 0; idx < cutPool->cnt; idx++) {
		duplicCut = FALSE;
		alpha = cutPool->vals[idx]->alpha;
		for (c = 0; c < cutsAdded->cnt; c++) {
			if (DBL_ABS(alpha - cutsAdded->vals[c]->alpha) < config.TOLERANCE) {
				if (equalVector(cutPool->vals[idx]->beta, cutsAdded->vals[c]->beta, betaLen, config.TOLERANCE)) {
					duplicCut = TRUE;
					break;
				}
			}
		}

		/* Add those cuts in cut pool that will be violated by incumbent solution */
		betaX = vXv(cutPool->vals[idx]->beta, incumbX, NULL, betaLen);
		if (betaX < alpha) {
			(*infeasIncumb) = TRUE;
			if (duplicCut == TRUE) {
				printf("Incumbent violates one old cut from feasible cut pool (this cut also exists in feasCutsAdded)\n");
				continue;
			}
			else
				printf( "Incumbent violates one new cut from feasible cut pool (this cut is not in feasCutsAdded but will be added)\n");
			cutsAdded->vals[cutsAdded->cnt++] = cutPool->vals[idx];

			printf("Cut added to master due to Incumbent violation\n");
		}
		else {
			/* Check if the cut will be violated by the candidate solution*/
			if (duplicCut == TRUE)
				continue;
			betaX = vXv(cutPool->vals[idx]->beta, candidX, NULL, betaLen);

			if (betaX < alpha) {
				printf("Candidate violates one cut from feasible cut pool (this cut is not in feasCutsAdded but will be added)\n");
				cutsAdded->vals[cutsAdded->cnt++] = cutPool->vals[idx];

				printf("Cut added to master due to candidate violation\n");
			}
		}
	}

	return 0;
}//END checkFeasCutPool()

/* This function will add a new feasibility cut to the master problem. Unlike addCut(), we do not rearrange the cuts while adding. */
int addfCut2Master(LPptr lp, oneCut *cut, vector incumbX, int lenX, int optCuts, int idx) {
	intvec 	indices;
	int		cnt;

	if (!(indices = (intvec) arr_alloc(lenX+1, int)))
		errMsg("Allocation", "addCut", "coefCol",0);
	for (cnt = 0; cnt < lenX; cnt++)
		indices[cnt + 1] = cnt;
	indices[0] = lenX;

	if ( config.MASTER_TYPE == PROB_QP )
		cut->alphaIncumb = cut->alpha - vXv(cut->beta, incumbX, NULL, lenX);

	/* add the row in the solver */
	if ( addRow(lp, lenX+1, cut->alphaIncumb, GE, 0, indices, cut->beta) ) {
		errMsg("solver", "addCut", "failed to add new row to problem in solver", 0);
		return 1;
	}
	cut->rowNum = lenX + optCuts + idx;

	mem_free(indices);
	return 0;
}//END addfCut()

void freeOneCut(oneCut *cut) {

	if (cut) {
		if (cut->iStar)
			mem_free(cut->iStar);
		if (cut->beta)
			mem_free(cut->beta);
		mem_free(cut);
	}
}

void freeCutsType(cutsType *cuts) {
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++)
		freeOneCut(cuts->vals[cnt]);
	mem_free(cuts->vals);
	mem_free(cuts);
}//END freeCuts
