/*
 * cuts.c
 *
 *  Created on: Apr 2, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "sdlp.h"

extern configType config;

int formCut(LPptr lp, LPptr sda, cellType *cell, probType *prob, cutsType *cuts, vector xt,
		int numRows, int numCols, BOOL isTerminal, int numStages, vector pi, intvec incumbCuts) {
	oneCut 	*cut;
	int		idxCut, status;

	/* allocate memory to new cut */
	cut = newCut(cell->omega->cnt, cell->k, prob->num->cntCcols);
	if ( cut == NULL ) {
		errMsg("algorithm", "formCut", "failed to allocate memory to the new cut", 0);
		return -1;
	}

	/* compute cut coefficients */
	status = stageCut(prob->num, prob->coord, cell->sigma, cell->delta, cell->omega, cell->lb, xt, cell->k, cut, isTerminal, numStages,
			cell->piRatios, &cell->dualStableFlag);
	if (status ) {
		errMsg("algorithm", "formCut", "failed to create the stage cut", 0);
		return -1;
	}

	/* add cut to cuts structure, decision simulation and stage dual approximation problems for previous stage */
	idxCut = addCut(lp, sda, cuts, numRows, numCols, cell->k, prob->num->cntCcols, prob->coord->colsC, cut,
			pi, incumbCuts);
	if ( idxCut < 0 ) {
		errMsg("algorithm", "formCut", "failed to add the cut stage problem", 0);
		freeOneCut(cut); return -1;
	}

	return idxCut;
}//END formCut()

int formIncumbCut(cellType *cell, probType *prob, LPptr lp, LPptr sda, cutsType *cuts, vector incumbU,
		int numRows, int numCols, BOOL isTerminal, int numStages, vector pi, intvec incumbCuts) {
	double 	mubBar;
	int 	extraRows = 0, idxSigma, idxCut;
	BOOL	newSigmaFlag;

	/* solve the subproblem with incumbent state as input */
	computeEndoRHS(prob->bBar, prob->Cbar, incumbU, cell->rhs);
	if ( computeExoRHS(cell->sda, NULL, prob->coord, prob->num, cell->omega->vals[cell->omega->idx],
			incumbU, cell->rhs) ){
		errMsg("allocation", "formIncumbCut", "failed to change the right-hand side with uncertainty and state information", 0);
		return -1;
	}
	if ( dualUpdates(cell->sda, cell->sp->name, prob->num->rows+extraRows, prob->num->cols, cell->pi, &mubBar)) {
		errMsg("algorithm", "formIncumbCut","failed to complete d%ual updates", 0);
		return -1;
	}

	/* update all the stochastic components, indicate that the updates with respect to new node have been completed TODO: future value */
	if (isTerminal )
		idxSigma = stocUpdate(config.MAX_ITER, prob->num, prob->coord, prob->Cbar, prob->bBar, cell->pi, mubBar, 0.0,
				cell->lambda, cell->sigma, &newSigmaFlag, cell->delta, cell->omega, cell->k);
	else
		idxSigma = stocUpdate(config.MAX_ITER, prob->num, prob->coord, prob->Cbar, prob->bBar, cell->pi, mubBar, cell->cuts->vals[cell->incumb->cutidx[0]]->alpha,
				cell->lambda, cell->sigma, &newSigmaFlag, cell->delta, cell->omega, cell->k);

#ifdef STOC_CHECK
	double obj;
	obj = cell->sigma->vals[idxSigma].pib - vXv(cell->sigma->vals[idxSigma].piC, incumbU, prob->coord->colsC, prob->num->cntCcols);
	obj += cell->delta->vals[cell->sigma->lambdaIdx[idxSigma]][cell->omega->idx].pib - vXv(cell->delta->vals[cell->sigma->lambdaIdx[idxSigma]][cell->omega->idx].piC,
			cell->omega->vals[cell->omega->idx], prob->coord->rvCols, prob->num->rvColCnt);
	printf("Objective function estimate at incumbent solution = %lf\n", obj);
#endif

	/* form new incumbent cut */
	idxCut = formCut(lp, sda, cell, prob, cuts, incumbU, numRows, numCols, isTerminal, numStages, pi, incumbCuts);
	if ( idxCut < 0 ) {
		errMsg("algorithm", "formIncumbCut", "failed to add the candidate cut", 0);
		return 1;
	}
	cuts->vals[idxCut]->isIncumb = TRUE;

	return idxCut;
}//END formIncumbCut()

/* subroutine to the allocate memory to oneCut structure and initialize its elements with default values */
oneCut *newCut(int numIstar, int numObs, int betaLen){
	oneCut *cut;

	if ( !(cut = (oneCut *) mem_malloc (sizeof(oneCut))) )
		errMsg("allocation", "newCut", "cut", 0);

	cut->numObs    = numObs;
	cut->numIstar  = numIstar;

	if ( !(cut->iStar = (intvec) arr_alloc(numIstar, int)) )
		errMsg("allocation", "newCut", "cut->iStar", 0);
	if ( !(cut->beta = (vector) arr_alloc(betaLen+1, double)) )
		errMsg("allocation", "newCut", "cut->beta", 0);
	cut->alpha 	 = 0.0;
	cut->beta[0] = 1.0;

	cut->isIncumb = FALSE;

	return cut;
}//END newCut()

/* This function allocates memory for a new cut structure.  This entails the structure itself, and the _val_ array of oneCut pointers
 * inside the structure.  The actual oneCut structures are allocated according to the numBeta parameter, via calls to new_cut(). */
cutsType *newCuts(int maxCuts) {
    cutsType *cuts;

    if (!(cuts = (cutsType *) mem_malloc (sizeof(cutsType))))
        errMsg("allocation", "newCuts", "cuts",0);
    if (!(cuts->vals = (oneCut **) arr_alloc (maxCuts, oneCut)))
        errMsg("allocation", "newCuts", "oneCuts",0);
    cuts->cnt = 0;
    cuts->maxCuts = maxCuts;

    return cuts;
}//END newCuts

int stageCut(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, omegaType *omega, double lb,
		vector xt, int numObs, oneCut *cut, BOOL isTerminal, int numStages, vector piRatios, BOOL *dualStableFlag) {
	vector 	pixC, beta;
	double	variance, estWindow, estAll, argmaxAll, argmaxWindow;
	int		cnt, i;
	iType	iStar;
	BOOL 	piEvalFlag = FALSE;

	if ( !(pixC = (vector) arr_alloc(sigma->cnt, double)) )
		errMsg("allocation", "stageCut", "pixC", 0);
	if ( !(beta = (vector) arr_alloc(num->prevCols+1, double)) )
		errMsg("allocation", "stageCut", "beta", 0);

    /* Determine if dual stability needs to be checked */
    if (numStages == 2 && numObs > config.PI_EVAL_START && !(numObs % config.PI_CYCLE))
        piEvalFlag = TRUE;

	/* calculate \bar{C}_t^\top \pi and \bar{C}_0^\top \pi, one at a time */
	for ( cnt = 0; cnt < sigma->cnt; cnt++ )
		for (i = 1; i <= num->cntCcols; i++)
			pixC[cnt] += sigma->vals[cnt].piC[i] * xt[coord->colsC[i]];

	estAll = 0.0; estWindow = 0.0;
	for (cnt = 0; cnt < omega->cnt; cnt++) {
		/* For each observation, find the Pi which maximizes height at X. */
		if ( piEvalFlag ) {
			iStar   	= computeIstar(num, coord, sigma, delta, pixC, xt, cnt, numObs, lb, isTerminal, &argmaxAll, FALSE);
			computeIstar(num, coord, sigma, delta, pixC, xt, cnt, numObs, lb, isTerminal, &argmaxWindow, TRUE);
			estAll 		+= argmaxAll*omega->weights[cnt];
			estWindow 	+= argmaxWindow*omega->weights[cnt];
		}
		else
			iStar = computeIstar(num, coord, sigma, delta, pixC, xt, cnt, numObs, lb, isTerminal, &argmaxAll, FALSE);

		/* identify the best stochastic element for all observations */
		cut->iStar[cnt] = iStar.sigma;

		/* Average using these pi's to calculate the cut coefficients. Here we are multiplying the coefficients by the number of times the observations
		 * are encountered: Intercept term first */
		cut->alpha += (sigma->vals[iStar.sigma].pib + delta->vals[iStar.delta][cnt].pib)* omega->weights[cnt];

		/* Slope term next */
		for (i = 1; i <= num->cntCcols; i++)
			beta[coord->colsC[i]] += sigma->vals[iStar.sigma].piC[i] * omega->weights[cnt];
		for (i = 1; i <= num->rvColCnt; i++)
			beta[coord->rvCols[i]] += delta->vals[iStar.delta][cnt].piC[i] * omega->weights[cnt];
	}

	if (piEvalFlag == TRUE) {
		piRatios[numObs % config.SCAN_LEN] = estWindow/ estAll;
		if (numObs - config.PI_EVAL_START >= config.SCAN_LEN)
			variance = calcVariance(piRatios, config.SCAN_LEN);
		else
			variance = calcVariance(piRatios, numObs);

		if ((DBL_ABS(variance) >= .000002) || (piRatios[numObs % config.SCAN_LEN] < 0.95) )
			*dualStableFlag = FALSE;
		else
			*dualStableFlag = TRUE;
	}

	/* Divide the coefficients by the number of observations to obtain the probability estimate */
	cut->alpha = cut->alpha/numObs;
	for (i = 1; i <= num->cntCcols; i++)
		cut->beta[i] = beta[coord->colsC[i]]/numObs;

	mem_free(pixC);
	mem_free(beta);

	return 0;
}//END stageCut()

iType computeIstar(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, vector pixC, vector xt, int cnt, int numObs, double lb, BOOL isTerminal,
		double *argmax, BOOL piEval) {
	iType 	iStar;
	int 	n, m, deltaIdx, window;
	double	arg, t_over_k;

	/* if piEval is TRUE then compute argmax using only the pi's generated in the window, otherwise, use all the pi's */
	if ( piEval )
		window = (int) (0.9*numObs);
	else
		window = numObs;

	(*argmax) = -DBL_MAX;

	for (n = 0; n < sigma->cnt; n++) {
		if ( sigma->ck[n] <= window ) {
			deltaIdx = sigma->lambdaIdx[n];

			/* Start with (\pi^\top \bar{b}) + (\pi^\top x \tilde{\omega}) - (\bar{C}_t^\top \pi)*x_t */
			arg = sigma->vals[n].pib + delta->vals[deltaIdx][cnt].pib - pixC[n];

			/* Subtract (\tilde{C}_t^\top \pi)*u_t */
			for (m = 1; m <= num->rvColCnt; m++)
				arg -= delta->vals[deltaIdx][cnt].piC[m] * xt[coord->rvCols[m]];

			/* Weigh the older dual solutions by the iteration count, this is done for all non-terminal stages */
			if ( !(isTerminal) ) {
				t_over_k = ((double) sigma->ck[n]/ (double) numObs);
				arg = arg * t_over_k + (1 - t_over_k)*lb;
			}

			if (arg > (*argmax)) {
				(*argmax) = arg;
				iStar.sigma = n;
				iStar.delta = deltaIdx;
			}
		}
	}

	return iStar;
}//END computeIstar()

int addCut(LPptr lp, LPptr sda, cutsType *cuts, int numRows, int numCols, int numObs, int betaLen, intvec betaIndices, oneCut *cut,
		vector pi, intvec incumbCuts) {
	intvec	indices;
	int		n;

	if ( !(indices = (intvec) arr_alloc(betaLen+1, int)) )
		errMsg("allocation", "addCut", "indices", 0);
	for ( n = 1; n <= betaLen; n++ )
		indices[n] = betaIndices[n]-1;
	indices[0] = numCols;

	/* make sure there is room to add a new cut */
	if ( cuts->cnt >= cuts->maxCuts) {
		if ( reduceCuts(lp, sda, cuts, numObs, incumbCuts, pi) ) {
			errMsg("algorithm", "addCut", "failed to reduce cuts", 0);
			mem_free(indices); return -1;
		}
	}

	if ( addRow(lp, betaLen+1, cut->alpha, 'G', 0, indices, cut->beta) ) {
		errMsg("solver", "addCut", "failed to add the new cut to solver problem", 0);
		return -1;
	}
	if ( sda != NULL ) {
		if ( cuts->cnt >= 1 )
			/* remove the previous timestage dual approximation before adding the new one */
			if ( removeRow(sda, numRows, numRows) ) {
				errMsg("solver", "addCut", "failed to add the new cut to solver problem", 0);
				return -1;
			}

		/* add current approximation */
		if ( addRow(sda, betaLen+1, cut->alpha, 'G', 0, indices, cut->beta) ) {
			errMsg("solver", "addCut", "failed to add the new cut to solver problem", 0);
			return -1;
		}
	}

	cut->rowNum = numRows + cuts->cnt;
	cuts->vals[cuts->cnt] = cut;

	mem_free(indices);
	return cuts->cnt++;
}//END addCut()

double maxCutHeight(cutsType *cuts, double lb, int iter, intvec Ccols, int betaLen, vector xt) {
	double Sm, ht;
	int cnt;

	Sm = cutHeight(cuts->vals[0], lb, iter, Ccols, betaLen, xt);
	for (cnt = 1; cnt < cuts->cnt; cnt++) {
		ht = cutHeight(cuts->vals[cnt], lb, iter, Ccols, betaLen, xt);
		if (Sm < ht)
			Sm = ht;
	}

	return Sm;
}//END maxCutHeight()

double cutHeight(oneCut *cut, double lb, int numObs, intvec Ccols, int betaLen, vector xt) {
	double height;
	double t_over_k = ((double) cut->numObs/ (double) numObs);

	/* A cut is calculated as Alpha - Beta x X */
	height = cut->alpha - vXv(cut->beta, xt, Ccols, betaLen);

	/* Weight cut based on number of observations used to form it */
	height *= t_over_k;

	/* account for non-trivial lower bound. */
	height += (1 - t_over_k) * lb;

	return height;
}//END cutHeight

/* This function will remove the oldest cut whose corresponding dual variable is zero (thus, a cut which was slack in last solution). */
int reduceCuts(LPptr lp, LPptr sda, cutsType *cuts, int numObs, intvec incumbCuts, vector pi) {
	int minObs, oldestCut, idx;

	if ( sda != NULL ) {
		errMsg("algorithm", "reduceCuts", "no scheme to reduce cuts for non-root stages", 0);
		return 1;
	}

	minObs 	  = numObs;
	oldestCut = cuts->cnt;

	/* identify the oldest loose cut */
    for (idx = 0; idx < cuts->cnt; idx++) {
		if ( idx == incumbCuts[0] )
			/* avoid dropping incumbent cut*/
			continue;

		if (cuts->vals[idx]->numObs < minObs && DBL_ABS(pi[cuts->vals[idx]->rowNum + 1]) <= config.TOLERANCE ) {
			minObs = cuts->vals[idx]->numObs;
			oldestCut = idx;
		}
	}

	/* if the oldest loose cut is the most recently added cut */
	if ( oldestCut == cuts->cnt ) {
		errMsg("algorithm", "reduceCuts", "failed to identify any cuts to drop", 0);
		return 1;
	}

	/* drop the selected cut and swap the last cut into its place */
	if ( dropCut(lp, cuts, oldestCut, incumbCuts) ) {
		errMsg("algorithm", "reduceCuts", "failed to drop a cut", 0);
		return 1;
	}

	return 0;
}//END reduceCuts()

/* This function removes a cut from both the cutType structure and the master problem constraint matrix.  In the cuts->val array, the last
 * cut is swapped into the place of the exiting cut.  In the constraint matrix, the row is deleted, and the row numbers of all constraints
 * below it are decremented. */
int dropCut(LPptr lp, cutsType *cuts, int cutIdx, intvec incumbCuts) {

	int idx, status, deletedRow;

	deletedRow = cuts->vals[cutIdx]->rowNum;
	/* Get rid of the indexed cut */
	status = removeRow(lp, deletedRow, deletedRow);
	if ( status ) {
		errMsg("solver", "dropCut", "failed to remove a row from master problem", 0);
		return 1;
	}
	freeOneCut(cuts->vals[cutIdx]);

	/* Update the surviving cuts */
	cuts->vals[cutIdx] = cuts->vals[--cuts->cnt];
	for (idx = 0; idx < cuts->cnt; idx++)
		if (cuts->vals[idx]->rowNum > deletedRow)
			--cuts->vals[idx]->rowNum;

    /* if the swapped cut happens to be the incumbent cut, then update its index */
    if ( incumbCuts[0] == cuts->cnt )
    	incumbCuts[0] = cutIdx;

	return 0;
}//END dropCut()

void freeCutsType(cutsType *cuts) {
	int n;

	if (cuts->vals) {
		for ( n = 0; n < cuts->cnt; n++ )
			if ( cuts->vals[n]) freeOneCut(cuts->vals[n]);
		mem_free(cuts->vals);
	}
	mem_free(cuts);

}//END freeCutsType()

void freeOneCut(oneCut *cut) {

	if (cut->beta) mem_free(cut->beta);
	if (cut->iStar) mem_free(cut->iStar);
	mem_free(cut);

}//END freeOneCut()

