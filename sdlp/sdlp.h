/*
 * sdlp.h
 *
 *  Created on: Apr 2, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *     Contact: harsha@smu.edu
 *
 */

#ifndef SDLP_H_
#define SDLP_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"

#undef CELL_SETUP

#define		TRIVIAL		0
#define		NONTRIVIAL	1

typedef struct {
	int			MAX_ITER;
	int			MIN_ITER;
	double  	TOLERANCE;
	long long 	RUN_SEED;
	int			EVAL_FLAG;
	long long 	EVAL_SEED;
	double		EVAL_ERROR;
	double		OPT_GAP;
}configType;

typedef struct {
	int		cnt;
	intvec	weights;
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

/* sdlp.c */
void parseCmdLine(string probName);
int readConfig(string inputDir);

/* algo.c */
int algo(oneProblem *orig, stocType *stoc, timeType *tim);
int forwardPass(probType **prob, cellType **cell, vector observ, int numStages);
int backwardPass(probType **prob, cellType **cell, int numStages);
void computeEndoRHS(sparseVector *bBar, sparseMatrix *Cbar, vector candidU, vector rhs);
int computeExoRHS(LPptr lp, coordType *coord, numType *num, vector observ, vector candidut, vector endoRHS);
void printAlgoDetails(int item);void cleanupAlgo(probType **prob, cellType **cell, int T);
void printAlgoDetails(int item);

/* setup.c */
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType ***cell);
cellType **newCell(stocType *stoc, probType **prob, vector lb, int T);
void freeCellType (probType **prob, cellType **cell, int T);

/* stocupdt.c */
int calcOmega(omegastuff *omegas, omegaType *omega, vector observ);
omegaType *newOmega(int t, stocType *stoc, int numObs);
lambdaType *newLambda(int numIter);
sigmaType *newSigma(int numIter, int numPi);
deltaType *newDelta(int numObs);
void freeLambdaType(lambdaType *lambda, BOOL all);
void freeSigmaType(sigmaType *sigma, BOOL all);
void freeDeltaType(deltaType *delta, int numObs, BOOL all);
void freeOmegaType(omegaType *omega);

/* cuts.c */
void freeCutsType(cutsType *cuts);
void freeOneCut(oneCut *cut);

/* optimal.c */
BOOL optimal(probType **prob, cellType **cell, int T);

#endif /* SDLP_H_ */
