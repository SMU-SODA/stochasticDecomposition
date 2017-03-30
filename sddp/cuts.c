/*
 * cuts.c
 *
 *  Created on: Dec 7, 2015
 *      Author: gjharsha
 */

#include <sddp.h>

extern configType config;

int formOptCut(probType *prob, cellType *cell, intvec iStar, LPptr lp, int numRows, int numCols, cutsType *cuts, vector U) {
	oneCut 	*cut;
	int		i, obs, idxCut;

	cut = newCut(cell->k, cell->omega->cnt, prob->num->cntCcols);

	for (obs = 0; obs < cell->omega->cnt; obs++ ) {
		/* intercept term */
		cut->alpha += (cell->sigma->vals[iStar[obs]].pib + cell->delta->vals[obs].pib)*cell->omega->probs[obs];

		/* Slope term next */
		for (i = 1; i <= prob->num->cntCcols; i++)
			cut->beta[prob->coord->colsC[i]] += cell->sigma->vals[iStar[obs]].piC[i] * cell->omega->probs[obs];
		for (i = 1; i <= prob->num->rvColCnt; i++)
			cut->beta[prob->coord->rvCols[i]] += cell->delta->vals[iStar[obs]].piC[i] * cell->omega->probs[obs];
	}
	cut->iStar = iStar;

	i = 0;
	while ( i < cuts->cnt ){
		if ( cuts->vals[i]->alpha - cut->alpha < config.TOLERANCE ) {
			if ( equalVector(cuts->vals[i]->beta, cut->beta, prob->num->cntCcols, config.TOLERANCE) ){
				freeOneCut(cut);
				return i;
			}
		}
		i++;
	}

	/* add cut to cuts structure, decision simulation and stage dual approximation problems for previous stage */
	idxCut = addCut(lp, cuts, numRows, numCols, prob->num->cntCcols, prob->coord->colsC, cut);
	if ( idxCut < 0 ) {
		errMsg("algorithm", "formCandidCut", "failed to add the cut stage problem", 0);
		return -1;
	}

	return cuts->cnt;
}//END formOptCut()

int addCut(LPptr lp, cutsType *cuts, int numRows, int numCols, int betaLen, intvec betaIndices, oneCut *cut) {
	intvec	indices;
	int		n;

	if ( !(indices = (intvec) arr_alloc(betaLen+1, int)) )
		errMsg("allocation", "addCut", "indices", 0);
	for ( n = 1; n <= betaLen; n++ )
		indices[n] = betaIndices[n]-1;
	indices[0] = numCols;

	if ( addRow(lp, betaLen+1, cut->alpha, 'G', 0, indices, cut->beta) )
		errMsg("solver", "addCut", "failed to add the new cut to solver problem", 0);

	cut->rowNum = numRows + cuts->cnt;
	cuts->vals[cuts->cnt] = cut;

	mem_free(indices);
	return cuts->cnt++;
}//END addCut()

int computeIstar(numType *num, coordType *coord, lambdaType *lambda, sigmaType *sigma, omegaType *omega, deltaType *delta, int obs, vector xt) {
	vector 	pixC;
	double	pib, arg, argmax;
	int 	m, n, iStar;

	argmax = -DBL_MAX; iStar = 0;
	if ( !(pixC = (vector) arr_alloc(num->rvColCnt+1, double)) )
		errMsg("allocation", "calcDeltaRow", "pixC", 0);

	for (n = 0; n < sigma->cnt; n++) {
		pib = vXv(lambda->vals[sigma->lambdaIdx[n]], omega->vals[obs], NULL, num->rvRowCnt);

		for ( m = num->rvRowCnt+1; m <= num->numRV; m++ )
			pixC[m-num->rvRowCnt] = lambda->vals[sigma->lambdaIdx[n]][m]*omega->vals[obs][m];

		/* Start with (\pi^\top \bar{b}) - (\bar{C}_t^\top \pi)*x_t */
		arg = sigma->vals[n].pib - vXv(sigma->vals->piC, xt, coord->colsC, num->cntCcols) +
				pib - vXv(pixC, xt, coord->rvCols, num->rvColCnt);

		if (arg > argmax) {
			argmax = arg;
			iStar = n;
		}
	}

	delta->vals[obs].pib = vXv(lambda->vals[sigma->lambdaIdx[iStar]], omega->vals[obs], NULL, num->rvRowCnt);
	for ( m = num->rvRowCnt+1; m <= num->numRV; m++ )
		pixC[num->rvRowCnt] = lambda->vals[sigma->lambdaIdx[iStar]][m]*omega->vals[obs][m];
	delta->vals[obs].piC = pixC;

	mem_free(pixC);
	return iStar;
}//END computeIstar()

/* subroutine to the allocate memory to oneCut structure and initialize its elements with default values */
oneCut *newCut(int ck, int numIstar, int betaLen) {
	oneCut *cut;

	if ( !(cut = (oneCut *) mem_malloc (sizeof(oneCut))) )
		errMsg("allocation", "newCut", "cut", 0);
	if ( !(cut->beta = (vector) arr_alloc(betaLen+1, double)) )
		errMsg("allocation", "newCut", "cut->beta", 0);

	cut->ck    = ck;
	cut->alpha 	 = 0.0;
	cut->beta[0] = 1.0;

	return cut;
}//END newCut()

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
