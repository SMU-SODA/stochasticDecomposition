/*
 * benders.h
 *
 *  Created on: Sep 21, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#ifndef BENDERS_H_
#define BENDERS_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"

#undef ALGO_CHECK
#undef STOCH_CHECK
#define SAVE_DUALS

typedef struct{
	int		NUM_REPS;			/* Maximum number of replications that can be carried out. */
	long long *RUN_SEED;		/* seed used during optimization */
	double 	TOLERANCE; 			/* for zero identity test */
	int		MIN_ITER;			/* minimum number of iterations */
	int		MAX_ITER;			/* maximum number of iterations */
	int		MASTER_TYPE;		/* type of master problem */
	double	MIN_QUAD_SCALAR;	/* Minimum value for regularizing parameter */
	double	EPSILON;			/* Optimality gap */

	int		EVAL_FLAG;
	int		NUM_EVALS;
	long long *EVAL_SEED;
	int		EVAL_MIN_ITER;
	double	EVAL_ERROR;

	int		CUT_MULT;			/* Determines the number of cuts to be used for approximate */
	double 	MAX_QUAD_SCALAR;	/* Maximum value for regularizing parameter */
	int		MULTICUT;			/* Set to 1 if multicut is to be solved */
	double	R1;
	double	R2;
	double	R3;

	int		MAX_OBS;			/* Maximum number of iterations before which SAA is invoked */
	int		SAA; 				/* Use SAA when continuous distribution in stoch file (1), or not (0) */

	int 	MULTIPLE_REP;		/* When multiple replications are needed, set this to (1), else (0) */
}configType;

typedef struct {
	int		ck;					/* Iteration when the cut was generated */
	double  alpha;              /* scalar value for the right-hand side */
	dVector  beta;               /* coefficients of the master problems's primal variables */
	bool	isIncumb;			/* indicates if the cut is an incumbent cut */
	double 	alphaIncumb;		/* right-hand side when using QP master, this is useful for quick updates */
	int 	rowNum;				/* row number for master problem in solver */
	cString	name;
}oneCut;

typedef struct {
	int    	cnt;                    /* number of cuts */
	oneCut  **vals;					/* values which define the set of cuts */
}cutsType;

typedef struct {
	double	repTime;
	double 	iterTime;
	double 	masterIter;
	double 	subprobIter;
	double 	optTestIter;
	double 	iterAccumTime;
	double 	masterAccumTime;
	double 	subprobAccumTime;
	double 	optTestAccumTime;
}runTime;

typedef struct {
	int		numRV;					/* Number of random variables */
	int 	cnt;					/* Number of observations */
	dVector	probs;					/* Probability of observation */
	dVector	*vals;					/* Observation values */
} omegaType;

typedef struct {
	int         k;                  /* number of iterations */
	int 		LPcnt; 				/* the number of LPs solved. */

    oneProblem  *master;            /* store master information */
	oneProblem 	*subprob;			/* store subproblem information */

	dVector      candidX;            /* primal solution of the master problem */
	double      candidEst;          /* objective value master problem */

	dVector      incumbX;			/* incumbent master solution */
	double      incumbEst;			/* estimate at incumbent solution */
	double 		quadScalar; 		/* the proximal parameter/quadratic scalar 'sigma' */
	bool        incumbChg;			/* set to be true if the incumbent solution has changed in an iteration */
	int         iCutIdx;			/* index of incumbent cut in cell->cuts structure */
	dVector		piM;

    int      	maxCuts;            /* maximum number of cuts to be used*/
	cutsType    *cuts;              /* optimality cuts */
	cutsType    *fCuts;             /* feasibility cuts */

	omegaType 	*omega;				/* all realizations observed during the algorithm */

    bool        optFlag;
    bool		optMode;
    bool		spFeasFlag;			/* Indicates whether the subproblem is feasible */
    int			feasCnt;			/* keeps track of the number of times infeasible candidate solution was encountered */
	bool		infeasIncumb;		/* indicates if the incumbent solution is infeasbible */

	runTime		*time;				/* Run time structure */
}cellType;

#if defined(SAVE_DUALS)
typedef struct {
	int 	cnt;
	dVector 	*vals;
	iVector	iter;
	iVector  obs;
}dualsType;
#endif

int parseCmdLine(int argc, char *argv[], cString probName, cString inputDir);
void createOutputDir(cString outputDir, cString algoName, cString probName);
int readConfig();
void freeConfig();

/* algo.c */
int algo (oneProblem *orig, timeType *tim, stocType *stoc, cString probName);
int solveBendersCell(stocType *stoc, probType **prob, cellType *cell);
bool optimal(cellType *cell);
void writeStatistic(FILE *soln, FILE *incumb, probType **prob, cellType *cell);

/* setup.c */
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell, dVector *meanSol);
cellType *newCell(stocType *stoc, probType **prob, dVector xk);
int cleanCellType(cellType *cell, probType *prob, dVector xk);
void freeCellType(cellType *cell);

/* masters.c */
int solveMaster(numType *num, sparseVector *dBar, cellType *cell);
int addCut2Master(cellType *cell, cutsType *cuts, oneCut *cut, int lenX);
int checkImprovement(probType *prob, cellType *cell, int candidCut);
int replaceIncumbent(probType *prob, cellType *cell);
int constructQP(probType *prob, cellType *cell, dVector incumbX, double quadScalar);
int changeQPproximal(LPptr lp, int numCols, double sigma);
int changeQPrhs(probType *prob, cellType *cell, dVector xk);
int changeQPbds(LPptr lp, int numCols, dVector bdl, dVector bdu, dVector xk);
oneProblem *newMaster(oneProblem *orig, double lb);

/* cuts.c */
int formOptCut(probType *prob, cellType *cell, dVector Xvect, bool isIncumb);
double maxCutHeight(cutsType *cuts, dVector xk, int betaLen);
double cutHeight(oneCut *cut, dVector xk, int betaLen);
int reduceCuts(oneProblem *master, cutsType *cuts, dVector vectX, dVector piM, int betaLen, int *iCutIdx);
int dropCut(oneProblem *master, cutsType *cuts, int cutIdx, int *iCutIdx);
oneCut *newCut(int numX, int currentIter);
cutsType *newCuts(int maxCuts);
void freeOneCut(oneCut *cut);
void freeCutsType(cutsType *cuts, bool partial);

/* subprob.c */
int solveSubprob(probType *prob, oneProblem *subproblem, dVector Xvect, dVector obsVals, bool *spFeasFlag, double *subprobTime, dVector piS, double *mubBar);
dVector computeRHS(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, dVector X, dVector obs);
dVector computeCostCoeff(numType *num, coordType *coord, sparseVector *dBar, dVector obs);
int computeMU(LPptr lp, int numCols, double *mubBar);
oneProblem *newSubproblem(oneProblem *subprob);
void chgRHSwSoln(sparseVector *bBar, sparseMatrix *Cbar, dVector rhs, dVector X);
int chgRHSwObserv(LPptr lp, numType *num, coordType *coord, dVector observ, dVector spRHS, dVector X);
int chgObjxwObserv(LPptr lp, numType *num, coordType *coord, dVector cost, iVector indices, dVector observ);
omegaType *newOmega(stocType *stoc);
void freeOmegaType(omegaType *omega, bool partial);

/* evaluate.c */
int evaluate(FILE *soln, stocType *stoc, probType **prob, cellType *cell, dVector Xvect);

#endif /* BENDERS_H_ */
