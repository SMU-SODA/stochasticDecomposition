/*
 * sddp.h
 *
 *  Created on: Dec 4, 2015
 *      Author: Harsha Gangammanavar
 */

#ifndef SDDP_H_
#define SDDP_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"

#undef CELL_SETUP
#undef STOC_CHECK
#undef CUT_CHECK
#define ALGO_RUN

#define		TRIVIAL		0
#define		NONTRIVIAL	1

typedef struct {
	int			MAX_ITER;
	int			MIN_ITER;
	double		TOLERANCE;
	long long	FORWPASS_SEED;
	double		BACKPASS_PCT;
	long long	BACKPASS_SEED;
	int			EVAL_FLAG;
	long long	EVAL_SEED;
	double		EVAL_ERROR;
	double		OPT_GAP;
}configType;

typedef struct {
	string	type;
	int		cnt;
	vector	probs;
	vector	*vals;
}omegaType;

typedef struct {
	int 	cnt;
	vector  *vals;
}lambdaType;

typedef struct {
	double	pib;
	vector  piC;
}pixbCType;

typedef struct {
	int 	cnt;
	intvec	lambdaIdx;
	pixbCType	*vals;
}sigmaType;

typedef struct {
	pixbCType	*vals;
}deltaType;

typedef struct {
	int delta;
	int sigma;
} iType;

typedef struct {
	int		ck;
	intvec	iStar;
	double	alpha;
	vector 	beta;
	int		rowNum;
}oneCut;

typedef struct {
	int		cnt;
	oneCut	**vals;
}cutsType;

/*
 * This is main structure which is used in the algorithm. This structure holds all information necessary to describe the problem at any particular
 * stage, in any iteration. Some elements of this structure are static/deterministic in nature, while others are dynamically updated (both in size as
 * well as value). The detailed description of each element is provided next to their declaration below.
 *
 */
typedef struct {
	int			k;				/* number of iterations */
	double		lb;				/* a lower bound computed using mean values for random variables */
	int			lbType;			/* TRIVIAL if lower bound is zero, and NONTRIVIAL if it is nonzero */
	oneProblem 	*sp;			/* stage subproblem, a copy of probType stage problem */
	vector		rhs;			/* right-hand side after state information update */
	vector		candidU;		/* candidate solution for the stage */
	vector		*incumbU;		/* a list of incumbent solutions for the stage. Used only when regularization is employed. */
	double		candidEst;		/* objective function value at candidate */
	vector		pi;				/* dual solution for the stage */
	cutsType	*cuts;			/* optimality cuts added */
	omegaType	*omega;			/* structure to hold observations */
	lambdaType	*lambda;		/* structure to hold dual information with for rows with random variables (roght-hand side or transfer matrix) */
	sigmaType	*sigma;			/* structure to hold "dual multiplied by deterministic part" */
	deltaType	*delta;			/* structure to hold "dual multiplied by stochastic part */
}cellType;

/* subroutines in sddp.c */
void parseCmdLine(string probName);
int readConfig(string inputDir);
void printAlgoDetails(int item);
void printSolutionDetails(int numStages, int t, probType **prob, cellType **cell);

/* subroutines in algo.c */
void printProbDetails(probType **p);
int algo (oneProblem *orig, stocType *stoc, timeType *tim);
int forwardPass(stocType *stoc, probType **prob, cellType **cell, int numStages);
int backwardPass(probType **prob, cellType **cell, int numStages);
void printProbDetails(probType **p);

/* subroutines in setup.c */
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType ***cell);
cellType **newCell(stocType *stoc, probType **prob, vector lb, int T);
vector calcLowerBound(oneProblem *orig, timeType *tim);
void cleanupAlgo(probType **prob, cellType **cell, int T);
void freeCellType (probType **prob, cellType **cell, int T);
void freeCutsType(cutsType *cuts);
void freeOneCut(oneCut *cut);

/* subroutine stocUpdt.c */
int stocUpdate(probType *prob, cellType *cell, vector candidU, int obs);
int computeMu(LPptr lp, int numCols, double *mubBar);
int calcLambda(numType *num, coordType *coord, lambdaType *lambda, vector pi, BOOL *newLambdaFlag);
int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, cutsType *cuts, vector pi, double mubBar, BOOL newLambdaFlag,
		int idxLambda, sigmaType *sigma);
void calcDelta(numType *num, vector lambdaPi, omegaType *omega, int obs, deltaType *delta);
omegaType *newOmega(int t, stocType *stoc);
lambdaType *newLambda(int numIter);
sigmaType *newSigma(int numIter, int numPi);
deltaType *newDelta(int numObs);
void freeLambdaType(lambdaType *lambda, BOOL all);
void freeSigmaType(sigmaType *sigma, BOOL all);
void freeDeltaType(deltaType *delta, int numObs, int all);
void freeOmegaType(omegaType *omega);

/* subroutines in modfiy.c */
void computeEndoRHS(sparseVector *bBar, sparseMatrix *Cbar, vector candidU, vector rhs);
int computeExoRHS(LPptr lp, coordType *coord, numType *num, vector observ, vector candidut, vector rhs);

/* subroutine in cuts.c */
oneCut *newCut(int ck, int numIstar, int betaLen);
int computeIstar(numType *num, coordType *coord, lambdaType *lambda, sigmaType *sigma, omegaType *omega, deltaType *delta, int obs, vector xt);
int formOptCut(probType *prob, cellType *cell, intvec iStar, LPptr lp, int numRows, int numCols, cutsType *cuts, vector U);
int addCut(LPptr lp, cutsType *cuts, int numRows, int numCols, int betaLen, intvec betaIndices, oneCut *cut);

/* subroutine in evaluate.c */
BOOL optimal(stocType *stoc, probType **prob, cellType **cell, int numStages);

#endif /* SDDP_H_ */
