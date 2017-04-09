/*
 * stocupdt.c
 *
 *  Created on: Dec 7, 2015
 *      Author: gjharsha
 */

#include <sddp.h>

extern configType config;

int stocUpdate(probType *prob, cellType *cell, vector candidU, int obs) {
	double 	mubBar;
	int		idxLambda, idxSigma;
	BOOL	newLambdaFlag;

	/* obtain the dual solution */
	if (getDual(cell->sp->lp, cell->pi, cell->sp->mar) ) {
		errMsg("solver", "backwardPass", "failed to obtain optimal dual solutions to stage problem", 0);
		return -1;
	}

	/* compute \bar{\mu} */
	if (computeMu(cell->sp->lp, cell->sp->mac, &mubBar) ) {
		errMsg("algorithm", "backwardPass", "failed to compute mu for stochastic updates", 0);
		return -1;
	}

#ifdef STOC_CHECK
	printf("Objective function value = %lf\t", getObjective(cell->sp->lp, PROB_LP));
#endif

	/* extract dual solutions corresponding to random rows into lambdaType */
	idxLambda = calcLambda(prob->num, prob->coord, cell->lambda, cell->pi, &newLambdaFlag);

	/* update the dual information with respect to bBar and Cbar */
	idxSigma = calcSigma(prob->num, prob->coord, prob->bBar, prob->Cbar,  cell->cuts, cell->pi, mubBar, newLambdaFlag, idxLambda, cell->sigma);

	calcDelta(prob->num, cell->lambda->vals[idxLambda], cell->omega, obs, cell->delta);

#ifdef STOC_CHECK
	double obj;
	obj = cell->sigma->vals[idxSigma].pib - vXv(cell->sigma->vals[idxSigma].piC, candidU, prob->coord->colsC, prob->num->cntCcols);
	obj += cell->delta->vals[obs].pib - vXv(cell->delta->vals[obs].piC, cell->omega->vals[obs], prob->coord->rvCols, prob->num->rvColCnt);
	printf("Objective function estimate = %lf\n", obj);
#endif

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

int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, cutsType *cuts, vector pi, double mubBar, BOOL newLambdaFlag,
		int idxLambda, sigmaType *sigma) {
	vector 	piCBar, temp;
	double 	pibBar;
	int 	cnt;

	/* sigma = (\pi_t^\top \bar{b}_t + \theta_t^\top \alpha_{t+}, \bar{C}_t^\top \pi_t) */
	pibBar = vXvSparse(pi, bBar) + mubBar;
	if ( cuts != NULL )
		for (cnt = 0; cnt < cuts->cnt; cnt++ )
			pibBar += cuts->vals[cnt]->alpha*pi[cuts->vals[cnt]->rowNum+1];

	temp = vxMSparse(pi, CBar, num->prevCols);
	piCBar = reduceVector(temp, coord->colsC, num->cntCcols);
	mem_free(temp);

	if (!newLambdaFlag){
		for (cnt = 0; cnt < sigma->cnt; cnt++) {
			if (DBL_ABS(pibBar - sigma->vals[cnt].pib) <= config.TOLERANCE) {
				if (equalVector(piCBar, sigma->vals[cnt].piC, num->cntCcols, config.TOLERANCE)) {
					mem_free(piCBar);
					return cnt;
				}
			}
		}
	}

	sigma->vals[sigma->cnt].pib = pibBar;
	sigma->vals[sigma->cnt].piC = piCBar;
	sigma->lambdaIdx[sigma->cnt] = idxLambda;

	return sigma->cnt++;
}//END calcSigma()

void calcDelta(numType *num, vector lambdaPi, omegaType *omega, int obs, deltaType *delta) {
	vector	pixC;
	int		m;

	if ( !(pixC = (vector) arr_alloc(num->rvColCnt+1, double)) )
		errMsg("allocation", "calcDeltaRow", "pixC", 0);

	delta->vals[obs].pib = 0;
	for ( m = 1; m <= num->rvRowCnt; m++ )
		delta->vals[obs].pib += lambdaPi[m]*omega->vals[obs][m];

	for ( m = num->rvRowCnt+1; m <= num->numRV; m++ )
		pixC[m-num->rvRowCnt] = lambdaPi[m]*omega->vals[obs][m];
	delta->vals[obs].piC = pixC;

}//END calcDeltaCol()

int computeMu(LPptr lp, int numCols, double *mubBar) {
	vector	dj, u;
	intvec	cstat;
	int		n;

	(*mubBar) = 0.0;

	if ( !(dj = (vector) arr_alloc(numCols+1, double)))
		errMsg("allocation", "computeMu", "dual slacks", 0);
	if ( !(u = (vector) arr_alloc(numCols+1, double)))
		errMsg("allocation", "computeMu", "TDA solutions", 0);

	if ( getPrimal(lp, u, numCols) ) {
		errMsg("solver", "forOptPass", "failed to obtain primal solution", 0);
		return 1;
	}
	if (getDualSlacks(lp, dj, numCols) ) {
		errMsg("solver", "computeMu", "failed to obtain dual slacks", 0);
		return 1;
	}

	/* extra column for eta if the stage problem is a QP */
	if ( !(cstat = (intvec) arr_alloc(numCols+2, int)) )
		errMsg("allocation", "computeMu", "column status", 0);
	if (getBasis(lp, cstat+1, NULL)) {
		errMsg("solver", "computeMu", "failed to get column status", 0);
		return 1;
	}

	for (n = 1; n <= numCols;  n++) {
		switch (cstat[n]) {
		case AT_LOWER:
			(*mubBar) += dj[n]*u[n];
			break;
		case AT_UPPER:
			(*mubBar) += dj[n]*u[n];
			break;
		default:
			break;
		}
	}

	mem_free(u); mem_free(cstat); mem_free(dj);

	return 0;
}//END computeMu()

omegaType *newOmega(int t, stocType *stoc) {
	omegaType *omega;
	int		  m, n;

	if ( !(omega = (omegaType *) mem_malloc(sizeof(omegaType))) )
		errMsg("allocation", "newOmega", "omega", 0);

	if ( strstr(stoc->type, "BLOCKS") ) {
		omega->cnt = stoc->numVals[t];
		if ( !(omega->probs = (vector) arr_alloc(omega->cnt, double)) )
			errMsg("allocation", "newOmega", "omega->probs", 0);
		if ( !(omega->vals = (vector *) arr_alloc(omega->cnt, vector)) )
			errMsg("allocation", "newOmega", "omega->vals", 0);

		for ( m = 0; m < omega->cnt; m++ ) {
			omega->probs[m] = stoc->probs[t][m];
			if ( !(omega->vals[m] = (vector) arr_alloc(stoc->numPerGroup[t]+1, double)) )
				errMsg("allocation", "newOmega", "omega->vals[m]", 0);
			for ( n = 0; n < stoc->numPerGroup[t]; n++ )
				omega->vals[m][n+1] = stoc->vals[stoc->groupBeg[t]+n][m] - stoc->mean[stoc->groupBeg[t]+n];
		}
	}
	else {
		errMsg("setup", "newOmega", "Non-Block structures are not currently supported", 0);
		return NULL;
	}


	return omega;
}//END newOmega()

lambdaType *newLambda(int numIter) {
	lambdaType *lambda;

	if ( !(lambda = (lambdaType *) mem_malloc(sizeof(lambdaType))) )
		errMsg("allocation", "newLambda", "cell->lambda", 0);
	if ( !(lambda->vals = (vector *) arr_alloc(numIter, vector)) )
		errMsg("allocation", "newLambda", "cell->lambda->vals", 0);
	lambda->cnt = 0;

	return lambda;
}//END newLambda()

sigmaType *newSigma(int numIter, int numPi) {
	sigmaType 	*sigma;

	if (!(sigma = (sigmaType *) mem_malloc(sizeof(sigmaType))) )
		errMsg("allocation", "newSigma", "sigma structure", 0);
	if ( !(sigma->lambdaIdx = (intvec) arr_alloc(numIter, int)) )
		errMsg("allocation", "newSigma", "sigma->lambdaIdx", 0);
	if ( !(sigma->vals = (pixbCType *) arr_alloc(numIter, pixbCType)) )
		errMsg("allocation", "newSigma", "sigma->vals", 0);
	sigma->cnt = numPi;

	return sigma;
}//END newSigma()

deltaType *newDelta(int numObs) {
	deltaType *delta;

	if ( !(delta = (deltaType *) mem_malloc(sizeof(deltaType))) )
		errMsg("allocation", "newDelta", "delta", 0);
	if ( !(delta->vals = (pixbCType *) arr_alloc(numObs, pixbCType)) )
		errMsg("allocation", "newDelta", "delta->vals", 0);

	return delta;
}//END newDelta()

void freeOmegaType(omegaType *omega) {
	int n;

	if (omega->probs) mem_free(omega->probs);
	if (omega->vals) {
		for ( n = 0; n < omega->cnt; n++ )
			if ( omega->vals[n] ) mem_free(omega->vals[n]);
		mem_free(omega->vals);
	}
	mem_free(omega);

}//END freeOmegaType()

void freeLambdaType(lambdaType *lambda, BOOL all) {
	int n;

	if (all) {
		mem_free(lambda->vals);
		mem_free(lambda);
	}
	else {
		for ( n = 0; n < lambda->cnt; n++ )
			if (lambda->vals[n]) mem_free(lambda->vals[n]);
		lambda->cnt = 0;
	}

}//END freeLambdaType()

void freeSigmaType(sigmaType *sigma, BOOL all) {
	int n;

	if ( all ) {
		if (sigma->vals) mem_free(sigma->vals);
		if (sigma->lambdaIdx) mem_free(sigma->lambdaIdx);
		mem_free(sigma);
	}
	else {
		for ( n = 0; n < sigma->cnt; n++ )
			if (sigma->vals[n].piC) mem_free(sigma->vals[n].piC);
		sigma->cnt = 0;
	}

}//END freeSigmaType()

void freeDeltaType(deltaType *delta, int numObs, int all) {
	int n;

	if (all) {
		mem_free(delta->vals);
		mem_free(delta);
	}
	else {
		for ( n = 0; n < numObs; n++ )
			if (delta->vals[n].piC) mem_free(delta->vals[n].piC);
	}

}//END freeDeltaType()
