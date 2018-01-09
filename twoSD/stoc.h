/*
 * stoc.h
 *
 *  Created on: Sep 11, 2017
 *      Author: gjharsha
 */

#ifndef STOC_H_
#define STOC_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"

#undef ALGO_CHECK
#undef STOCH_CHECK

/* To save time and space, Pi x b and Pi x C are calculated as soon as possible and stored in structures like sigma and delta.  Towards this end,
 * pixbCType represents a single calculation of pi X b (which is a scalar) and pi X C (which is a vector).*/
typedef struct{
	double 	pib;					/* scalar pi x b */
	vector 	piC;					/* vector pi x C */
} pixbCType;

/* The OmegaType data structure stores the set of observations which have been made so far. Each observation consists of a vector
 * of realizations _vals_ of random variables. Each rv occurs in right-hand side vector, cost coefficient vector or the T matrix
 * which (along with the candidate first-stage solution) make up the rhs of the subproblem. The _weight_ field specifies the
 * number of times a particular outcome has been observed (it starts at 1 when the outcome is first generated, and increments
 * every time the same outcome is observed again).  _cnt_ just specifies the number of distinct outcomes which have been observed
 ** and stored in the omegaType structure. */
typedef struct {
	int		numRV;
	int 	cnt;
	intvec	weights;                 /* number of times that an omega is observed */
	vector	probs;
	vector	*vals;
} omegaType;

/* The lambda structure stores some of the dual variable values from every distinct dual vector obtained during the program.  Each vector contains
 * only those dual variables whose corresponding rows in the subproblem constraint matrix contain random elements.  _val_ is an array of
 * these dual vectors (thus it is 2-D). _row_ gives the corresponding row number for a given dual variable in _val_.  _cnt_ represents
 * the number of dual vectors currently stored in lambda. */
typedef struct {
	int 	cnt;
	vector 	*vals;
} lambdaType;

/* The sigma matrix contains the values of Pi x bBar and Pi x Cbar  for all values of pi obtained so far (note it does not depend on
 * observations of omega).  _col_ gives the column number of each non-zero element in pi X Cbar.  _val_ is an array of values
 * for pi X bBar and pi X Cbar, one entry for each pi.  Note that values  which are always zero (because Rbar or Cbar is zero there) are not
 * stored.  The _lamb_ array is the same size as the _val_ array, and for each element in _val_ the corresponding element in _lamb_ references
 * the dual vector in lambda that was used to calculate that entry in sigma. */
typedef struct {
	int 		cnt;
	pixbCType 	*vals;
	intvec		lambdaIdx;
	intvec		ck;
} sigmaType;

/* The delta matrix contains the values of lambda_pi X bOmega and lambda_pi X Comega for all values of pi and all observations of omega.
 * _col_ gives the column number of the each non-zero element in the multiplication of lambda_pi X Comega (the same elements are non-zero
 * each time).  _val_ is an array of vectors of (lambda_pi X bOmega, lambda_pi X Comega) pairs of calculations.  A row in _val_ corresponds
 * to a distinct dual vector, and a column in _val_ corresponds to a distinct observation of omega.  Thus, every pi-omega combination is
 * represented here, and the size of the delta matrix can be determine from lambda->cnt and omega->cnt.
 * Note that when elements of omega get dropped, vacant columns appear in delta.  This is ok, but be sure to loop carefully! */
typedef struct {
	pixbCType 	**vals;
} deltaType;

/* subprob.c */
int solveSubprob(probType *prob, oneProblem *subproblem, vector Xvect, lambdaType *lambda, sigmaType *sigma, deltaType *delta, int deltaRowLength,
		omegaType *omega, int omegaIdx, BOOL *newOmegaFlag, int currentIter, double TOLERANCE, BOOL *subFeasFlag, BOOL *newSigmaFlag,
		double *subprobTime, double *argmaxTime);
vector computeRHS(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, vector X, vector obs);
void chgRHSwSoln(sparseVector *bBar, sparseMatrix *Cbar, vector rhs, vector X) ;
int chgRHSwObserv(LPptr lp, numType *num, coordType *coord, vector observ, vector spRHS, vector X);
oneProblem *newSubprob(oneProblem *sp);

/* stocUpdate.c */
int stochasticUpdates(probType *prob, LPptr spLP, lambdaType *lambda, sigmaType *sigma, deltaType *delta, int deltaRowLength, omegaType *omega,
		int omegaIdx, BOOL newOmegaFlag, int currentIter, double TOLERANCE);
int computeIstar(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, vector piCbarX, vector Xvect, int obs,
		int numSamples, BOOL pi_eval, double *argmax, BOOL isNew);
void calcDelta(numType *num, coordType *coord, lambdaType *lambda, omegaType *omega, deltaType *delta, int deltaRowLength, int elemIdx,
		BOOL newOmegaFlag);
int calcLambda(numType *num, coordType *coord, vector Pi, lambdaType *lambda, BOOL *newLambdaFlag, double TOLERANCE);
int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, vector pi, double mubBar,
              int idxLambda, BOOL newLambdaFlag, int currentIter, sigmaType *sigma, BOOL *newSigmaFlag, double TOLERANCE);
int calcOmega(vector observ, int begin, int end, omegaType *omega, BOOL *newOmegaFlag, double TOLERANCE);
int computeMU(LPptr lp, intvec cstat, int numCols, double *mubBar);
lambdaType *newLambda(int num_iter, int numLambda, int numRVrows);
sigmaType *newSigma(int numIter, int numNzCols, int numPi);
deltaType *newDelta(int numIter);
omegaType *newOmega(int numIter);
void freeLambdaType(lambdaType *lambda, BOOL partial);
void freeSigmaType(sigmaType *sigma, BOOL partial);
void freeOmegaType(omegaType *omega, BOOL partial);
void freeDeltaType (deltaType *delta, int lambdaCnt, int omegaCnt, BOOL partial);

#endif /* STOC_H_ */
