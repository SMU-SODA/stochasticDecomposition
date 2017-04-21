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
#define ALGO_RUN
#undef STOC_CHECK

#define		TRIVIAL		0
#define		NONTRIVIAL	1

typedef struct {
	int			MAX_ITER;
	int			MIN_ITER;
	double  	TOLERANCE;
	long long 	RUN_SEED;
	int			QUADRATIC;
	double		MIN_QUAD_SCALAR;
	double		MAX_QUAD_SCALAR;
	double		R1;
	double		R2;
	double		R3;
	int			POLICY;
	int			EVAL_FLAG;
	long long 	EVAL_SEED;
	double		EVAL_ERROR;
	double		OPT_GAP;
}configType;

typedef struct {
	int		cnt;
	int		idx;
	BOOL	newObs;
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
	int 		cnt;
	intvec		ck;
	intvec		lambdaIdx;
	pixbCType	*vals;
}sigmaType;

typedef struct {
	pixbCType	**vals;
}deltaType;

/* When calculating istar for a cut, it is useful to have two separate references into the sigma and delta structures, since each dual vector
 * is stored in two places -- part in sigma and part in delta.  The final entry in cut->istar[] will just be the _sigma_ field of this structure. */
typedef struct {
	int delta;
	int sigma;
} iType;

typedef struct {
	int		numObs;
	int		numIstar;
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
	int		cnt;
	vector	*vals;
	vector	est;
	double  quadScalar;
	double	improv;
	double	normd_k;
	double  normd_k_1;
	BOOL	chg;
}incumbType;

typedef struct {
	int			k;				/* number of iterations */
	double		lb;				/* a lower bound computed using mean values for random variables */
	int			lbType;			/* TRIVIAL if lower bound is zero, and NONTRIVIAL if it is nonzero */
	oneProblem 	*sp;			/* stage subproblem used for decision simulation, a copy of probType stage problem */
	void		*sda;			/* stage dual approximation problem pointer to be used by solver */
	vector		rhs;			/* right-hand side after state information update */
	vector		candidU;		/* candidate solution for the stage */
	double		candidEst;		/* objective function estimate at candidate */
	incumbType	*incumb;		/* structure to hold all relevant information regarding incumbent solutions */
	vector		pi;				/* dual solution for the stage */
	int			maxCuts;		/* maximum number of cuts to be included in the cost-to-go function approximation */
	cutsType	*cuts;			/* optimality cuts added */
	omegaType	*omega;			/* structure to hold observations */
	lambdaType	*lambda;		/* structure to hold dual information with for rows with random variables (right-hand side or transfer matrix) */
	sigmaType	*sigma;			/* structure to hold "dual multiplied by deterministic part" */
	deltaType	*delta;			/* structure to hold "dual multiplied by stochastic part */
}cellType;

/* sdlp.c */
void parseCmdLine(string probName);
int readConfig(string inputDir);

/* algo.c */
int algo(string probName, oneProblem *orig, stocType *stoc, timeType *tim);
int forwardPass(probType **prob, cellType **cell, vector observ, int numStages);
int backwardPass(probType **prob, cellType **cell, vector observ, int numStages);
void computeEndoRHS(sparseVector *bBar, sparseMatrix *Cbar, vector candidU, vector rhs);
int computeExoRHS(LPptr lp, LPptr sda, coordType *coord, numType *num, vector observ, vector candidut, vector rhs);
int changeEtaCol(LPptr lp, int numCols, int numRows, int k, cutsType *cuts, double lb);
int updateCutsRHS(LPptr lp, cutsType *cuts, double lb, int numObs);
int dualUpdates(LPptr lp, string name, int numRows, int numCols, vector pi, double *mubBar);
int computeMu(LPptr lp, int numCols, double *mubBar);
void printAlgoDetails(int item);
void printSolutionDetails (string probName, probType **prob, cellType **cell, int numStages);
void cleanupAlgo(probType **prob, cellType **cell, int T);

/* setup.c */
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType ***cell);
cellType **newCell(stocType *stoc, probType **prob, vector lb, vector meanSol, int T);
incumbType *newIncumb(int maxIncumb, vector meanSol, int lenX, sparseVector *dBar);
void freeCellType (probType **prob, cellType **cell, int T);

/* stocupdt.c */
int stocUpdate(int maxIter, numType *num, coordType *coord, sparseMatrix *Cbar, sparseVector *bBar, vector pi, double mubBar, double futureVal,
		lambdaType *lambda, sigmaType *sigma, BOOL *newSigmaFlag, deltaType *delta, omegaType *omega, int numObs);
int calcOmega(omegastuff *omegas, omegaType *omega, vector observ);
int calcLambda(numType *num, coordType *coord, lambdaType *lambda, vector pi, BOOL *newLambdaFlag);
int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, vector pi, double mubBar, double futureVal,
		int idxLambda, BOOL newLambdaFlag, int numObs, sigmaType *sigma, BOOL *newSigmaFlag);
void calcDeltaCol(numType *num, coordType *coord, lambdaType *lambda, omegaType *omega, deltaType *delta);
void calcDeltaRow(int numIter, numType *num, coordType *coord, lambdaType *lambda, int idxLambda, omegaType *omega, deltaType *delta);
omegaType *newOmega(int t, stocType *stoc, int numObs);
lambdaType *newLambda(int numIter);
sigmaType *newSigma(int numIter, int numPi);
deltaType *newDelta(int numObs);
void freeLambdaType(lambdaType *lambda);
void freeSigmaType(sigmaType *sigma);
void freeDeltaType(deltaType *delta, int numObs, int numLambda);
void freeOmegaType(omegaType *omega);

/* cuts.c */
int formCandidCut(LPptr lp, LPptr sda, cellType *cell, probType *prob, cutsType *cuts, vector xt,
		int numRows, int numCols, int maxCuts, BOOL isTerminal);
oneCut *newCut(int numIstar, int numObs, int betaLen);
int stageCut(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, omegaType *omega,
		vector xt, int numObs, oneCut *cut, BOOL isTerminal);
int addCut(LPptr lp, LPptr sda, cutsType *cuts, int numRows, int numCols, int maxCuts, int betaLen, intvec betaIndices, oneCut *cut);
iType computeIstar(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, vector pixC, vector xt, int cnt, int numObs, BOOL isTerminal);
double maxCutHeight(cutsType *cuts, double lb, int iter, intvec Ccols, int betaLen, vector xt);
double cutHeight(oneCut *cut, double lb, int numObs, intvec Ccols, int betaLen, vector xt);
void freeCutsType(cutsType *cuts);
void freeOneCut(oneCut *cut);

/* quad.c */
int constructQP(LPptr lp, int numCols, double regSigma);
int changeQPrhs(LPptr lp, intvec betaCols, int betaLen, int numRows, sparseMatrix *Dbar, sparseVector *bBar, cutsType *cuts, vector X, vector rhs,
		int numObs, double lb);
int changeQPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector X);

/* policy.c */
int selectIncumb(incumbType *incumb);
void checkImprovement(probType **prob, cellType **cell, int numStages);

/* optimal.c */
BOOL optimal(probType **prob, cellType **cell, int T);

#endif /* SDLP_H_ */
