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

	for ( n = 0; n < omegas->numRV; n++ )
		observ[n] -= omegas->mean[n];

	n = 0;
	while ( n < omega->cnt ) {
		if ( equalVector(observ, omega->vals[n], omegas->numRV, config.TOLERANCE) )
			break;
		n++;
	}
	if ( n == omega->cnt ) {
		/* new observation encountered, store its values */
		if ( !(omega->vals[n] = (vector) arr_alloc(omegas->numRV+1,double)) )
			errMsg("allocation", "forwardPass", "cell[t]->omega->vals[n]", 0);
		copyVector(observ, omega->vals[n], omegas->numRV, 0);
		omega->weights[omega->cnt] = 1;
		return omega->cnt++;
	}

	omega->weights[n]++;

	return n;
}//END calcOmega

omegaType *newOmega(int t, stocType *stoc, int numObs) {
	omegaType *omega = NULL;

	if ( !(omega = (omegaType *) mem_malloc(sizeof(omegaType))) )
		errMsg("allocation", "newOmega", "omega", 0);
	if ( !(omega->vals = (vector *) arr_alloc(numObs, vector)) )
		errMsg("allocation", "newOmega", "omega->vals", 0);
	if ( !(omega->weights = (intvec) arr_alloc(numObs, double)) )
		errMsg("allocation", "newOmega", "omega->probs", 0);
	omega->cnt = 0;

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
	if ( !(sigma->vals = (pixbCType *) arr_alloc(numIter, pixbCType)) )
		errMsg("allocation", "newSigma", "sigma->vals", 0);
	sigma->cnt = numPi;

	return sigma;
}//END newSigma()

deltaType *newDelta(int numObs) {
	deltaType *delta = NULL;

	if ( !(delta = (deltaType *) mem_malloc(sizeof(deltaType))) )
		errMsg("allocation", "newDelta", "delta", 0);
	if ( !(delta->vals = (pixbCType *) arr_alloc(numObs, pixbCType)) )
		errMsg("allocation", "newDelta", "delta->vals", 0);

	return delta;
}//END newDelta()

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

void freeDeltaType(deltaType *delta, int numObs, BOOL all) {
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
