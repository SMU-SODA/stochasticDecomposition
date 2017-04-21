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

int formCandidCut(LPptr lp, LPptr sda, cellType *cell, probType *prob, cutsType *cuts, vector xt,
		int numRows, int numCols, int maxCuts, BOOL isTerminal) {
	oneCut 	*cut;
	int		idxCut, status;

	/* allocate memory to new cut */
	cut = newCut(cell->omega->cnt, cell->k, prob->num->cntCcols);
	if ( cut == NULL ) {
		errMsg("algorithm", "formCandidCut", "failed to allocate memory to the new cut", 0);
		return -1;
	}

	/* compute cut coefficients */
	status = stageCut(prob->num, prob->coord, cell->sigma, cell->delta, cell->omega, xt, cell->k, cut, isTerminal);
	if (status ) {
		errMsg("algorithm", "formNewCut", "failed to create the stage cut", 0);
		return -1;
	}

	/* add cut to cuts structure, decision simulation and stage dual approximation problems for previous stage */
	idxCut = addCut(lp, sda, cuts, numRows, numCols, maxCuts, prob->num->cntCcols, prob->coord->colsC, cut);
	if ( idxCut < 0 ) {
		errMsg("algorithm", "formCandidCut", "failed to add the cut stage problem", 0);
		return -1;
	}

	return idxCut;
}//END formCandidCut()

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

    return cuts;
}//END newCuts

int stageCut(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, omegaType *omega,
		vector xt, int numObs, oneCut *cut, BOOL isTerminal) {
	vector 	pixC, beta;
	int		cnt, i;
	iType	iStar;

	if ( !(pixC = (vector) arr_alloc(sigma->cnt, double)) )
		errMsg("allocation", "stageCut", "pixC", 0);
	if ( !(beta = (vector) arr_alloc(num->prevCols, double)) )
		errMsg("allocation", "stageCut", "beta", 0);

	/* calculate \bar{C}_t^\top \pi and \bar{C}_0^\top \pi, one at a time */
	for ( cnt = 0; cnt < sigma->cnt; cnt++ )
		for (i = 1; i <= num->cntCcols; i++)
			pixC[cnt] += sigma->vals[cnt].piC[i] * xt[coord->colsC[i]];

	for (cnt = 0; cnt < omega->cnt; cnt++) {
		/* For each observation, find the Pi which maximizes height at X. */
		iStar = computeIstar(num, coord, sigma, delta, pixC, xt, cnt, numObs, isTerminal);

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

	/* Divide the coefficients by the number of observations to obtain the probability estimate */
	cut->alpha = cut->alpha/numObs;
	for (i = 1; i <= num->cntCcols; i++)
		cut->beta[i] = beta[coord->colsC[i]]/numObs;

	mem_free(pixC);
	mem_free(beta);

	return 0;
}//END stageCut()

iType computeIstar(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, vector pixC, vector xt, int cnt, int numObs, BOOL isTerminal) {
	iType 	iStar;
	int 	n, m, deltaIdx;
	double	arg, argmax;

	argmax = -DBL_MAX;

	for (n = 0; n < sigma->cnt; n++) {
		deltaIdx = sigma->lambdaIdx[n];

		/* Start with (\pi^\top \bar{b}) + (\pi^\top x \tilde{\omega}) - (\bar{C}_t^\top \pi)*x_t */
		arg = sigma->vals[n].pib + delta->vals[deltaIdx][cnt].pib - pixC[n];

		/* Subtract (\tilde{C}_t^\top \pi)*u_t */
		for (m = 1; m <= num->rvColCnt; m++)
			arg -= delta->vals[deltaIdx][cnt].piC[m] * xt[coord->rvCols[m]];

		/* Weigh the older dual solutions by the iteration count, this is done for all non-terminal stages */
		if ( !(isTerminal) )
			arg = (arg * ((double) sigma->ck[n]))/(double) numObs;

		if (arg > argmax) {
			argmax = arg;
			iStar.sigma = n;
			iStar.delta = deltaIdx;
		}
	}

	return iStar;
}//END computeIstar()

int addCut(LPptr lp, LPptr sda, cutsType *cuts, int numRows, int numCols, int maxCuts, int betaLen, intvec betaIndices, oneCut *cut) {
	intvec	indices;
	int		n;

	if ( !(indices = (intvec) arr_alloc(betaLen+1, int)) )
		errMsg("allocation", "addCut", "indices", 0);
	for ( n = 1; n <= betaLen; n++ )
		indices[n] = betaIndices[n]-1;
	indices[0] = numCols;

	/* make sure there is room to add a new cut */
	if (cuts->cnt >= maxCuts) {
		errMsg("algorithm", "addCut", "ran out of memory for cuts", 0);
		return -1;
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

