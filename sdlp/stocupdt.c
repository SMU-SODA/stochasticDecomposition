/*
 * sddp.c
 *
 *  Created on: Mar 27, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *     Contact: harsha@smu.edu
 *
 */

#include "sdlp.h"

extern configType config;

int calcOmega(omegastuff *omegas, omegaType *omega, vector observ) {
	int n;

	for ( n = 1; n <= omegas->numRV; n++ )
		observ[n-1] -= omegas->mean[n];

	n = 0;
	while ( n < omega->cnt ) {
		if ( equalVector(observ-1, omega->vals[n], omegas->numRV, config.TOLERANCE) )
			break;
		n++;
	}
	if ( n == omega->cnt ) {
		/* new observation encountered, store its values */
		if ( !(omega->vals[n] = (vector) arr_alloc(omegas->numRV+1,double)) )
			errMsg("allocation", "forwardPass", "cell[t]->omega->vals[n]", 0);
		copyVector(observ, omega->vals[n], omegas->numRV, FALSE);
		omega->vals[n][0] = oneNorm(omega->vals[n]+1, omegas->numRV);
		omega->weights[omega->cnt] = 1;
		omega->newObs = TRUE;
		return omega->cnt++;
	}

	omega->weights[n]++;
	omega->newObs = FALSE;

	return n;
}//END calcOmega

int stocUpdate(int maxIter, numType *num, coordType *coord, sparseMatrix *Cbar, sparseVector *bBar, vector pi, double mubBar, double futureVal,
		lambdaType *lambda, sigmaType *sigma, BOOL *newSigmaFlag, deltaType *delta, omegaType *omega, int numObs) {
	int 	idxLambda, idxSigma;
	BOOL	newLambdaFlag;

	/* Only need to calculate column if new observation of omega found */
	if (omega->newObs)
		calcDeltaCol(num, coord, lambda, omega, delta);

	/* extract dual solutions corresponding to random rows into lambdaType */
	idxLambda = calcLambda(num, coord, lambda, pi, &newLambdaFlag);

	/* update the dual information with respect to bBar and Cbar */
	idxSigma = calcSigma(num, coord, bBar, Cbar, pi, mubBar, futureVal, idxLambda, newLambdaFlag, numObs, sigma, newSigmaFlag);

	/* need to calculate new row only if new lambda is observed in lambdaType */
	if (newLambdaFlag)
		calcDeltaRow(maxIter, num, coord, lambda, idxLambda, omega, delta);

	return idxSigma;
}//END stocUpdate()

int calcLambda(numType *num, coordType *coord, lambdaType *lambda, vector pi, BOOL *newLambdaFlag) {
	vector 	newPi;
	int		n, len;

	len = num->rvRowCnt;

	if (!(newPi = (vector) arr_alloc(len+1, double)) )
		errMsg("allocation", "calcLambda", "current dual vector corresponding to random rows", 0);

	for ( n = 1; n <= len; n++ )
		newPi[n] = pi[coord->rvRows[n]];

	/* compare new dual vector with previously observed duals stored in lambda structure */
	for ( n = 0; n < lambda->cnt; n++)
		if (equalVector(newPi, lambda->vals[n], len, config.TOLERANCE)) {
			mem_free(newPi);
			*newLambdaFlag = FALSE;
			return n;
		}

	/* add the vector to lambda structure */
	lambda->vals[lambda->cnt] = newPi;
	*newLambdaFlag = TRUE;

	return lambda->cnt++;
}//END calcLambda()

int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, vector pi, double mubBar, double futureVal,
		int idxLambda, BOOL newLambdaFlag, int numObs, sigmaType *sigma, BOOL *newSigmaFlag) {
	double pibBar;
	vector	piCBar, temp;
	int		cnt;

	/* sigma = \pi_t^\top \bar{b}_t - \bar{C}_t^\top \pi_t */
	pibBar = vXvSparse(pi, bBar) + mubBar + futureVal;

	temp = vxMSparse(pi, CBar, num->prevCols);
	piCBar = reduceVector(temp, coord->colsC, num->cntCcols);
	mem_free(temp);

	if (!newLambdaFlag){
		for (cnt = 0; cnt < sigma->cnt; cnt++) {
			if (DBL_ABS(pibBar - sigma->vals[cnt].pib) <= config.TOLERANCE) {
				if (equalVector(piCBar, sigma->vals[cnt].piC, num->cntCcols, config.TOLERANCE)) {
					mem_free(piCBar);
					(*newSigmaFlag) = FALSE;
					return cnt;
				}
			}
		}
	}

	(*newSigmaFlag) = TRUE;
	sigma->vals[sigma->cnt].pib  = pibBar;
	sigma->vals[sigma->cnt].piC  = piCBar;
	sigma->lambdaIdx[sigma->cnt] = idxLambda;
	sigma->ck[sigma->cnt] = numObs;

	return sigma->cnt++;
}//END calcSigma()

void calcDeltaCol(numType *num, coordType *coord, lambdaType *lambda, omegaType *omega, deltaType *delta) {
	int		idxLambda, m, idxOmega;
	vector	pixC;

	idxOmega = omega->idx;

	for ( idxLambda = 0; idxLambda < lambda->cnt; idxLambda++ ) {
		if ( num->rvColCnt > 0 ) {
			if ( !(pixC = (vector) arr_alloc(num->rvColCnt+1, double)) )
				errMsg("allocation", "calcDeltaCol", "pixC", 0);
		}
		else
			pixC = NULL;

		delta->vals[idxLambda][idxOmega].pib = 0;
		for ( m = 1; m <= num->rvRowCnt; m++ )
			delta->vals[idxLambda][idxOmega].pib += lambda->vals[idxLambda][m]*omega->vals[idxOmega][m];

		for ( m = num->rvRowCnt+1; m <= num->numRV; m++ )
			pixC[m-num->rvRowCnt] = lambda->vals[idxLambda][m]*omega->vals[idxOmega][m];
		delta->vals[idxLambda][idxOmega].piC = pixC;
	}

}//END calcDeltaCol()

void calcDeltaRow(int numIter, numType *num, coordType *coord, lambdaType *lambda, int idxLambda, omegaType *omega, deltaType *delta) {
	int 	idxOmega, m;
	vector 	pixC;
#ifdef TRACE
	trPrint("calcDeltaRow", 1);
#endif


	/* allocate memory to a new row in delta structure */
	if ( !(delta->vals[idxLambda] = (pixbCType *) arr_alloc(numIter, pixbCType)) )
		errMsg("allocation", "calcDeltaRow", "delta->vals[obs]", 0);

	/* For all observations, calculate \pi \times b and \pi \times C */
	for (idxOmega = 0; idxOmega < omega->cnt; idxOmega++) {
		if ( num->rvColCnt > 0 ) {
			if ( !(pixC = (vector) arr_alloc(num->rvColCnt+1, double)) )
				errMsg("allocation", "calcDeltaCol", "pixC", 0);
		}
		else
			pixC = NULL;

		delta->vals[idxLambda][idxOmega].pib = 0;
		for ( m = 1; m <= num->rvRowCnt; m++ )
			delta->vals[idxLambda][idxOmega].pib += lambda->vals[idxLambda][m]*omega->vals[idxOmega][m];

		for ( m = num->rvRowCnt+1; m <= num->numRV; m++ )
			pixC[m-num->rvRowCnt] = lambda->vals[idxLambda][m]*omega->vals[idxOmega][m];
		delta->vals[idxLambda][idxOmega].piC = pixC;
	}

	mem_free(pixC);
}//END calcDeltaRow()

omegaType *newOmega(int t, stocType *stoc, int numObs) {
	omegaType *omega = NULL;

	if ( !(omega = (omegaType *) mem_malloc(sizeof(omegaType))) )
		errMsg("allocation", "newOmega", "omega", 0);
	if ( !(omega->vals = (vector *) arr_alloc(numObs, vector)) )
		errMsg("allocation", "newOmega", "omega->vals", 0);
	if ( !(omega->weights = (intvec) arr_alloc(numObs, double)) )
		errMsg("allocation", "newOmega", "omega->probs", 0);
	omega->cnt = 0;
	omega->idx = 0;
	omega->newObs = FALSE;

	return omega;
}//END newOmega()

lambdaType *newLambda(int numIter) {
	lambdaType *lambda = NULL;

	if ( !(lambda = (lambdaType *) mem_malloc(sizeof(lambdaType))) )
		errMsg("allocation", "newLambda", "cell->lambda", 0);
	if ( !(lambda->vals = (vector *) arr_alloc(numIter, vector)) )
		errMsg("allocation", "newLambda", "cell->lambda->vals", 0);
	lambda->cnt = 0;

	return lambda;
}//END newLambda()

sigmaType *newSigma(int numIter, int numPi) {
	sigmaType 	*sigma = NULL;

	if (!(sigma = (sigmaType *) mem_malloc(sizeof(sigmaType))) )
		errMsg("allocation", "newSigma", "sigma structure", 0);
	if ( !(sigma->lambdaIdx = (intvec) arr_alloc(numIter, int)) )
		errMsg("allocation", "newSigma", "sigma->lambdaIdx", 0);
	if ( !(sigma->ck = (intvec) arr_alloc(numIter, int)) )
		errMsg("allocation", "newSigma", "sigma->ck", 0);
	if ( !(sigma->vals = (pixbCType *) arr_alloc(numIter, pixbCType)) )
		errMsg("allocation", "newSigma", "sigma->vals", 0);
	sigma->cnt = numPi;

	return sigma;
}//END newSigma()

deltaType *newDelta(int numIter) {
	deltaType *delta = NULL;

	if ( !(delta = (deltaType *) mem_malloc(sizeof(deltaType))) )
		errMsg("allocation", "newDelta", "delta", 0);
	if ( !(delta->vals = (pixbCType **) arr_alloc(numIter, pixbCType)) )
		errMsg("allocation", "newDelta", "delta->vals", 0);

	return delta;
}//END newDelta()

void freeLambdaType(lambdaType *lambda) {
	int n;

	if ( lambda->vals ) {
		for ( n = 0; n < lambda->cnt; n++ )
			if (lambda->vals[n]) mem_free(lambda->vals[n]);
		mem_free(lambda->vals);
	}
	mem_free(lambda);

}//END freeLambdaType()

void freeSigmaType(sigmaType *sigma) {
	int n;

	if (sigma->lambdaIdx) mem_free(sigma->lambdaIdx);
	if (sigma->ck) mem_free(sigma->ck);
	if ( sigma->vals) {
		for ( n = 0; n < sigma->cnt; n++ )
			if (sigma->vals[n].piC) mem_free(sigma->vals[n].piC);
		mem_free(sigma->vals);
	}
	mem_free(sigma);

}//END freeSigmaType()

void freeDeltaType(deltaType *delta, int numObs, int numLambda) {
	int m, n;

	if ( delta->vals ) {
		for ( n = 0; n < numLambda; n++ ) {
			if ( delta->vals[n] ) {
				for ( m = 0; m < numObs; m++ )
					if ( delta->vals[n][m].piC ) {
						mem_free(delta->vals[n][m].piC);
						delta->vals[n][m].piC = NULL;
					}
				mem_free(delta->vals[n]);
			}
		}
			mem_free(delta->vals);
	}
	mem_free(delta);

}//END freeDeltaType()

void freeOmegaType(omegaType *omega) {
	int n;

	if (omega->vals) {
		for (n = 0; n < omega->cnt; n++ )
			if (omega->vals[n]) mem_free(omega->vals[n]);
		mem_free(omega->vals);
	}
	if (omega->weights) mem_free(omega->weights);
	mem_free(omega);

}//END freeOmegaType
