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
#undef BASIS_CHECK

/* To save time and space, Pi x b and Pi x C are calculated as soon as possible and stored in structures like sigma and delta.  Towards this end,
 * pixbCType represents a single calculation of pi X b (which is a scalar) and pi X C (which is a dVector).*/
typedef struct{
	double 	pib;					/* scalar pi x b */
	dVector 	piC;					/* dVector pi x C */
} pixbCType;

/* The OmegaType data structure stores the set of observations which have been made so far. Each observation consists of a dVector
 * of realizations _vals_ of random variables. Each rv occurs in right-hand side dVector, cost coefficient dVector or the T matrix
 * which (along with the candidate first-stage solution) make up the rhs of the subproblem. The _weight_ field specifies the
 * number of times a particular outcome has been observed (it starts at 1 when the outcome is first generated, and increments
 * every time the same outcome is observed again).  _cnt_ just specifies the number of distinct outcomes which have been observed
 ** and stored in the omegaType structure. */
typedef struct {
	int		numRV;
	int 	cnt;
	iVector	weights;                 /* number of times that an omega is observed */
	dVector	probs;
	dVector	*vals;
} omegaType;

/* The lambda structure stores some of the dual variable values from every distinct dual dVector obtained during the program.  Each dVector contains
 * only those dual variables whose corresponding rows in the subproblem constraint matrix contain random elements.  _val_ is an array of
 * these dual dVectors (thus it is 2-D). _row_ gives the corresponding row number for a given dual variable in _val_.  _cnt_ represents
 * the number of dual dVectors currently stored in lambda. */
typedef struct {
	int 	cnt;
	dVector 	*vals;
} lambdaType;

/* The sigma matrix contains the values of Pi x bBar and Pi x Cbar  for all values of pi obtained so far (note it does not depend on
 * observations of omega).  _col_ gives the column number of each non-zero element in pi X Cbar.  _val_ is an array of values
 * for pi X bBar and pi X Cbar, one entry for each pi.  Note that values  which are always zero (because Rbar or Cbar is zero there) are not
 * stored.  The _lamb_ array is the same size as the _val_ array, and for each element in _val_ the corresponding element in _lamb_ references
 * the dual dVector in lambda that was used to calculate that entry in sigma. */
typedef struct {
	int 		cnt;
	pixbCType 	*vals;
	iVector		lambdaIdx;
	iVector		ck;
} sigmaType;

/* The delta matrix contains the values of lambda_pi X bOmega and lambda_pi X Comega for all values of pi and all observations of omega.
 * _col_ gives the column number of the each non-zero element in the multiplication of lambda_pi X Comega (the same elements are non-zero
 * each time).  _val_ is an array of dVectors of (lambda_pi X bOmega, lambda_pi X Comega) pairs of calculations.  A row in _val_ corresponds
 * to a distinct dual dVector, and a column in _val_ corresponds to a distinct observation of omega.  Thus, every pi-omega combination is
 * represented here, and the size of the delta matrix can be determine from lambda->cnt and omega->cnt.
 * Note that when elements of omega get dropped, vacant columns appear in delta.  This is ok, but be sure to loop carefully! */
typedef struct {
	pixbCType 	**vals;
} deltaType;

/* This structure is used to hold information about the sample used in the current iteration. Subproblems corresponding to the observations
 * in this sample are solved to optimality, and the solutions are used to conduct the stochastic updates.
 */
typedef struct {
	int		cnt;
	iVector	omegaIdx;				/* Observation index in omegaType */
	bool 	*newOmegaFlag;  		/* Flag indicates if the observation is encountered for the first time. */
	iVector basisIdx;	    		/* Basis index in basisType */
	bool	*newBasisFlag;  		/* Flag indicates if the basis is encountered for the first time. */
}sampleType;

typedef struct {
	int				ck;			/* The first time the basis was encountered. */
	int				weight;		/* Frequency of observation for each unique basis */
	unsigned long 	*rCode;		/* Encoded row status in the basis (currently not being used */
	unsigned long	*cCode;		/* Encoded column status in the basis */
	int				phiLength;	/* Number of basic columns with random cost coefficients */
	dVector			*phi;		/* The phi matrix: the columns of inverse dual basis matrix which have random cost coefficients */
	iVector			omegaIdx;	/* Indices within the random cost coefficient dVector to which the columns of phi matrix correspond to. */
	iVector			sigmaIdx;	/* Indices within the random cost coefficient dVector to which the columns of phi matrix correspond to. */
	dVector			piDet;		/* Deterministic component of the dual solution. This depends only on the basis. */
	double			mubBar;
	dVector			gBar;
	sparseMatrix	*psi;		/* The simplex tableau matrix corresponding to the basis. */
	bool			feasFlag;
}oneBasis;

/* The basis type data structure holds all the information regarding the basis identified during the course of the algorithm.
 * This structure will be at the heart of all calculations related to stochastic updates. */
typedef struct {
	int			basisDim;	/* The dimension of the basis matrix */
	int			cnt;		/* Number of unique basis encountered by the algorithm */
	int         init;       /* first index for evaluating the basis structure in the argmax */
	int			rCodeLen;	/* Length of encoded row status */
	int			cCodeLen;	/* Length of encoded column status */
	bool		**obsFeasible;
	oneBasis	**vals;		/* a structure for each basis */
}basisType;

/* subprob.c */
int solveSubprob(probType *prob, oneProblem *subproblem, dVector Xvect, basisType *basis, lambdaType *lambda, sigmaType *sigma, deltaType *delta, int deltaRowLength,
		omegaType *omega, int omegaIdx, bool *newOmegaFlag, int currentIter, double TOLERANCE, bool *subFeasFlag, bool *newBasisFlag,
		double *subprobTime, double *argmaxTime);
int computeRHS(LPptr lp, numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, dVector X, dVector obs);
int computeCostCoeff(LPptr lp, numType *num, coordType *coord, sparseVector *dBar, dVector observ);
void chgRHSwSoln(sparseVector *bBar, sparseMatrix *Cbar, dVector rhs, dVector X) ;
int chgRHSwObserv(LPptr lp, numType *num, coordType *coord, dVector observ, dVector spRHS, dVector X);
int chgObjxwObserv(LPptr lp, numType *num, coordType *coord, dVector cost, iVector indices, dVector observ);
oneProblem *newSubprob(oneProblem *sp);

/* stocUpdate.c */
int stochasticUpdates(probType *prob, LPptr spLP, basisType *basis, lambdaType *lambda, sigmaType *sigma, deltaType *delta, int deltaRowLength,
		omegaType *omega, int omegaIdx, bool newOmegaFlag, int currentIter, double TOLERANCE, bool *newBasisFlag, bool subFeasFlag);
int computeIstar(numType *num, coordType *coord, basisType *basis, sigmaType *sigma, deltaType *delta, sampleType *sample,
		dVector piCbarX, dVector Xvect, dVector observ, int obs, int numSamples, bool pi_eval, double *argmax, bool isNew);
int calcDelta(numType *num, coordType *coord, lambdaType *lambda, deltaType *delta, int deltaRowLength, omegaType *omega, bool newOmegaFlag, int elemIdx);
int calcLambda(numType *num, coordType *coord, dVector Pi, lambdaType *lambda, bool *newLambdaFlag, double TOLERANCE);
int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, dVector pi, double mubBar,
              int idxLambda, bool newLambdaFlag, int currentIter, sigmaType *sigma, bool *newSigmaFlag, double TOLERANCE);
int calcOmega(dVector observ, int begin, int end, omegaType *omega, bool *newOmegaFlag, double TOLERANCE);
int computeMU(LPptr lp, iVector cstat, int numCols, double *mubBar);
lambdaType *newLambda(int num_iter, int numLambda, int numRVrows);
sigmaType *newSigma(int numIter, int numNzCols, int numPi);
deltaType *newDelta(int numIter);
omegaType *newOmega(int numOmega, int numIter);
void freeLambdaType(lambdaType *lambda, bool partial);
void freeSigmaType(sigmaType *sigma, bool partial);
void freeOmegaType(omegaType *omega, bool partial);
void freeDeltaType (deltaType *delta, int lambdaCnt, int omegaCnt, bool partial);

/* randCost.c */
void calcBasis(LPptr lp, numType *num, coordType *coord, sparseVector *dBar, oneBasis *B, int basisDim);
int decomposeDualSolution(LPptr spLP, oneBasis *B, dVector omegaVals, int numRows);
oneBasis *newBasis(LPptr lp, int numCols, int numRows, int currentIter, bool subFeasFlag);
bool checkBasisFeasibility(oneBasis *B, sparseVector dOmega, cString senx, int numCols, int numRows, double TOLERANCE);
basisType *newBasisType(int numIter, int numCols, int numRows, int wordLength);
void freeOneBasis(oneBasis *B);
void freeBasisType(basisType *basis, bool partial);
sampleType *newSample(int sampleSize);
void freeSampleType(sampleType *sample);

#endif /* STOC_H_ */
