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

#include "stoc.h"
#include "twoSD.h"

extern configType config;

int formSDCut(probType *prob, cellType *cell, vector Xvect, int omegaIdx, BOOL newOmegaFlag, BOOL isIncumb) {
	oneCut *cut;
	int    cutIdx;

	/* (a) Construct the subproblem with input observation and master solution, solve the subproblem, and complete stochastic updates */
	if ( solveSubprob(prob, cell, Xvect, omegaIdx, newOmegaFlag) ) {
		errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
		return -1;
	}

	/* (b) create an affine lower bound */
	cut = SDCut(prob->num, prob->coord, cell->basis, cell->sigma, cell->delta, cell->omega, Xvect, cell->k, &cell->dualStableFlag, cell->pi_ratio, cell->lb);
	if ( cut == NULL ) {
		errMsg("algorithm", "formSDCut", "failed to create the affine minorant", 0);
		return -1;
	}

	/* (c) add cut to the master problem  */
	if ( (cutIdx = addCut2Master(cell, cut, prob->num->prevCols, cell->lb)) < 0 ) {
		errMsg("algorithm", "formSDCut", "failed to add the new cut to master problem", 0);
		return -1;
	}
	return cutIdx;
}//END formCut()

oneCut *SDCut(numType *num, coordType *coord, basisType *basis, sigmaType *sigma, deltaType *delta, omegaType *omega, vector Xvect, int numSamples,
		BOOL *dualStableFlag, vector pi_ratio, double lb) {
	oneCut *cut;
	vector 	piCbarX, beta;
	double  argmaxAll, argmaxNew, argmax, alpha = 0.0, argmax_dif_sum = 0.0, argmax_all_sum = 0.0, vari = 1.0;
	int 	istarAll, istarNew, istar, c, cnt, obs, offset;
	BOOL    pi_eval_flag = FALSE;

	/* allocate memory to hold a new cut */
	cut = newCut(num->prevCols, omega->cnt, numSamples);

	/* Need to store  Pi x Cbar x X independently of observation loop */
	if (!(piCbarX= arr_alloc(sigma->cnt, double)))
		errMsg("Allocation", "SDCut", "pi_Tbar_x",0);
	if ( !(beta = (vector) arr_alloc(num->prevCols + 1, double)) )
		errMsg("Allocation", "SDCut", "beta", 0);

	/* Calculate pi_eval_flag to determine the way of computing argmax */
	if (numSamples > config.PI_EVAL_START && !(numSamples % config.PI_CYCLE))
		pi_eval_flag = TRUE;

	/* Calculate (Pi x Cbar) x X by mult. each VxT by X, one at a time */
	for (cnt = 0; cnt < sigma->cnt; cnt++)
		piCbarX[cnt] = vXv(sigma->vals[cnt].piC, Xvect, coord->colsC, num->cntCcols);

	offset = num->rvbOmCnt + num->rvCOmCnt;
	/* Test for omega issues */
	for (obs = 0; obs < omega->cnt; obs++) {
		/* For each observation, find the Pi which maximizes height at X. */
		if (pi_eval_flag == TRUE) {
			istarAll = computeIstar(num, coord, basis, sigma, delta, Xvect, piCbarX, omega->vals[obs]+offset, obs, numSamples, pi_eval_flag, &argmaxAll, FALSE);
			istarNew = computeIstar(num, coord, basis, sigma, delta, Xvect, piCbarX, omega->vals[obs]+offset, obs, numSamples, TRUE, &argmaxNew, TRUE);

			if (argmaxNew > argmaxAll) {
				argmax = argmaxNew; istar = istarNew;
			}
			else {
				argmax = argmaxAll; istar = istarAll;
			}

			argmax_dif_sum += max(argmaxAll - lb, 0) * omega->weight[obs];
			argmax_all_sum += max(argmax - lb, 0) * omega->weight[obs];
		}
		else {
			/* identify the maximal Pi for each observation */
			istar = computeIstar(num, coord, basis, sigma, delta, Xvect, piCbarX, omega->vals[obs], obs, numSamples, pi_eval_flag, &argmax, FALSE);
		}

		if ( istar < 0 ) {
			errMsg("algorithm", "SDCut", "failed to identify maximal Pi for an observation", 0);
			return NULL;
		}
		cut->iStar[obs] = istar;

		/* Average using these Pi's to calculate the cut itself (update alpha and beta) */
		alpha += (sigma->vals[basis->vals[istar]->sigmaIdx[0]].pib + delta->vals[basis->vals[istar]->lambdaIdx[0]][obs].pib)* omega->weight[obs];

		for (c = 1; c <= num->cntCcols; c++)
			beta[coord->colsC[c]] += sigma->vals[basis->vals[istar]->sigmaIdx[0]].piC[c] * omega->weight[obs];
		for (c = 1; c <= num->rvCOmCnt; c++)
			beta[coord->rvCols[c]] += delta->vals[basis->vals[istar]->lambdaIdx[0]][obs].piC[c] * omega->weight[obs];

		for ( cnt = 1; cnt <= basis->vals[istar]->phiLength; cnt++ ) {
			alpha += (sigma->vals[basis->vals[istar]->sigmaIdx[cnt]].pib + delta->vals[basis->vals[istar]->lambdaIdx[cnt]][obs].pib)
					* (omega->vals[obs][offset+basis->vals[istar]->omegaIdx[cnt]]) * omega->weight[obs];
			for (c = 1; c <= num->cntCcols; c++)
				beta[coord->colsC[c]] += sigma->vals[basis->vals[istar]->sigmaIdx[cnt]].piC[c] * (omega->vals[obs][offset+basis->vals[istar]->omegaIdx[cnt]]) * omega->weight[obs];
			for (c = 1; c <= num->rvCOmCnt; c++)
				beta[coord->rvCols[c]] += delta->vals[basis->vals[istar]->lambdaIdx[cnt]][obs].piC[c] * (omega->vals[obs][offset+basis->vals[istar]->omegaIdx[cnt]]) * omega->weight[obs];
		}
	}

	if (pi_eval_flag == TRUE) {
		pi_ratio[numSamples % config.SCAN_LEN] = argmax_dif_sum / argmax_all_sum;
		if (numSamples - config.PI_EVAL_START > config.SCAN_LEN)
			vari = calcVari(pi_ratio, NULL, NULL, 0);

		if (DBL_ABS(vari) >= .000002 || (pi_ratio[numSamples % config.SCAN_LEN]) < 0.95)
			*dualStableFlag = FALSE;
		else
			*dualStableFlag = TRUE;
	}

	cut->alpha = alpha / numSamples;

	for (c = 1; c <= num->prevCols; c++)
		cut->beta[c] = beta[c] / numSamples;

	/* coefficient of eta coloumn */
	cut->beta[0] = 1.0;

	mem_free(piCbarX);
	mem_free(beta);

	return cut;
}//END SDCut

/*This function loops through all the dual vectors found so far and returns the index of the one which satisfies the expression:
 * 				argmax { Pi x (R - T x X) | all Pi }
 * where X, R, and T are given.  It is calculated in this form:
 * 				Pi x bBar + Pi x bomega + (Pi x Cbar) x X + (Pi x Comega) x X.
 * Since the Pi's are stored in two different structures (sigma and delta), the index to the maximizing Pi is actually a structure
 * containing two indices.  (While both indices point to pieces of the dual vectors, sigma and delta may not be in sync with one
 * another due to elimination of non-distinct or redundant vectors. */
int computeIstar(numType *num, coordType *coord, basisType *basis, sigmaType *sigma, deltaType *delta, vector Xvect, vector PiCbarX, vector omegaVals, int obs,
		int numSamples, BOOL pi_eval, double *argmax, BOOL isNew) {
	double 	arg;
	int 	sigmaIdx, lambdaIdx, cnt, i, maxCnt, basisUp, basisLow;

	if (pi_eval == TRUE)
		numSamples -= (numSamples / 10 + 1);

	/* Establish the range of iterations over which the istar calculations are conducted. Only bases discovered in this iteration range are used. */
	if ( !isNew ) {
		basisUp = numSamples; basisLow = -INT_MAX;
	}
	else {
		basisUp = INT_MAX; basisLow = numSamples;
	}

	*argmax = -DBL_MAX; maxCnt = 0;
	/* Run through the list of basis to choose the one which provides the best lower bound */
	for ( cnt = 0; cnt < basis->cnt; cnt++ ) {
		if ( basis->obsFeasible[cnt][obs] ) {
			/* I. Compute argument using deterministic component of the dual solution */
			sigmaIdx  = basis->vals[cnt]->sigmaIdx[0];
			lambdaIdx = basis->vals[cnt]->lambdaIdx[0];

			if ( basis->vals[cnt]->ck > basisLow && basis->vals[cnt]->ck <= basisUp ) {
				/* a. */
				arg = sigma->vals[sigmaIdx].pib + delta->vals[lambdaIdx][obs].pib - PiCbarX[sigmaIdx];

				/* b. */
				arg -= vXv(delta->vals[lambdaIdx][obs].piC, Xvect, coord->rvCols, num->rvCOmCnt);

				/* II. Append argument using values computed from the stochastic component of the dual solution. */
				for ( i = 1; i <= basis->vals[cnt]->phiLength; i++ ) {
					sigmaIdx  = basis->vals[cnt]->sigmaIdx[i];
					lambdaIdx = basis->vals[cnt]->lambdaIdx[i];

					/* a. */
					arg += omegaVals[basis->vals[cnt]->omegaIdx[i]]*(sigma->vals[sigmaIdx].pib + delta->vals[lambdaIdx][obs].pib - PiCbarX[sigmaIdx]);

					/* b. */
					arg -= omegaVals[basis->vals[cnt]->omegaIdx[i]]*vXv(delta->vals[lambdaIdx][obs].piC, Xvect, coord->rvCols, num->rvCOmCnt);
				}

				if (arg > (*argmax)) {
					*argmax = arg;
					maxCnt = cnt;
				}
			}
		}
	}

	return maxCnt;
}//END computeIstar

/* This function allocates memory for the arrays inside a single cut, and initializes its values accordingly.  The cut structure
 * itself is assumed to be already allocated.  Note, each beta vector contains room for its one-norm, thought it just gets filled
 * with zero anyway. */
oneCut *newCut(int numX, int numIstar, int numSamples) {
	oneCut *cut;

	cut = (oneCut *) mem_malloc (sizeof(oneCut));
	cut->cutObs   = numSamples;
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
int reduceCuts(cellType *cell, vector candidX, vector pi, int betaLen, double lb) {
	double height, minHeight;
	int minObs, oldestCut,idx;

	minObs 	  = cell->k;
	oldestCut = cell->cuts->cnt;

	/* identify the oldest loose cut */
	for (idx = 0; idx < cell->cuts->cnt; idx++) {
		if ( idx == cell->iCutIdx || cell->cuts->vals[idx]->rowNum < 0)
			/* avoid dropping incumbent cut and newly added cuts */
			continue;

		if (cell->cuts->vals[idx]->cutObs < minObs && DBL_ABS(pi[cell->cuts->vals[idx]->rowNum + 1]) <= config.TOLERANCE ) {
			minObs = cell->cuts->vals[idx]->cutObs;
			oldestCut = idx;
		}
	}

	/* if the oldest loose cut is the most recently added cut, then the cut with minimium cut height will be dropped */
	if ( oldestCut == cell->cuts->cnt ) {
		//minHeight = cutHeight(lbType, cell[agentIdx]->cuts->vals[0], cell[agentIdx]->k, candidX, betaLen, lb);
		minHeight = cutHeight(cell->cuts->vals[0], cell->k, candidX, betaLen, lb);
		oldestCut = 0;

		for (idx = 1; idx < cell->cuts->cnt; idx++) {
			if (idx == cell->iCutIdx)
				continue;

			//height = cutHeight(lbType, cell[agentIdx]->cuts->vals[idx], cell[agentIdx]->k, candidX, betaLen, lb);
			height = cutHeight(cell->cuts->vals[idx], cell->k, candidX, betaLen, lb);
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

	for (idx = 0; idx < cell->cuts->cnt; idx++) {
		if (cell->cuts->vals[idx]->rowNum > deletedRow)
			--cell->cuts->vals[idx]->rowNum;
	}

	/* decrease the number of rows on solver */
	cell->master->mar--;

	return 0;
}//END dropCut()

/*
 ** This function calculate the variance of the
 ** vector x.
 */
double calcVari(double *x, double *mean_value, double *stdev_value, int batch_size) {
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
