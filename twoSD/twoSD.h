/*
 * twoSD.h
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#ifndef TWOSD_H_
#define TWOSD_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"
#include "stoc.h"

#define TRIVIAL 0
#define NONTRIVIAL 1
#define INF	DBL_MAX

#undef STOCH_CHECK
#undef ALGO_CHECK

/* A data structure which holds on the configuration information about the algorithm. Most of these configuration parameters are read from a
configuration file. These elements, once set during initialization, are not modified during the course of the algorithm. */
typedef struct{
	long long RUN_SEED;			/* seed used during optimization */
	double 	TOLERANCE; 			/* for zero identity test */
	int		MIN_ITER;			/* minimum number of iterations */
	int		MAX_ITER;			/* maximum number of iterations */
	int		MASTERTYPE;			/* type of master problem */
	int		CUT_MULT;			/* Determines the number of cuts to be used for approximate */
	int		TAU;				/* Frequency at which the incumbent is updated */
	double	MIN_QUAD_SCALAR;	/* Minimum value for regularizing parameter */
	double 	MAX_QUAD_SCALAR;	/* Maximum value for regularizing parameter */
	double	R1;
	double	R2;
	double	R3;
	int		PI_EVAL_START;
	int		PI_CYCLE;
	int		SCAN_LEN;
	int		EVAL_FLAG;
	long long EVAL_SEED;
	int		EVAL_MIN_ITER;
	double	EVAL_ERROR;
	double  PRE_EPSILON;
	double	EPSILON;
	int		BOOTSTRAP_REP;		/* Number of boot-strap replications in full optimality test */
	double	PERCENT_PASS;		/* percentage of bootstrap replications need to be satisfied */
}configType;

/* The oneCut and cutsType data structures will be used to hold all information which can completely define the affine minorants (cuts) which
 * are used to compute the lower bounding function approximations */
typedef struct {
	double  alpha;                  /* scalar value for the right-hand side */
	vector  beta;                   /* coefficients of the master problems's primal variables */
	int 	cutObs;					/* number of samples on which the given cut was based */
	int 	omegaCnt;				/* number of *distinct* observations on which the cut is based (this is also the length of istar) */
	intvec	iStar;					/* indices of maximal pi for each distinct observation */
	BOOL	isIncumb;				/* indicates if the cut is an incumbent cut */
	double 	alphaIncumb;			/* right-hand side when using QP master, this is useful for quick updates */
	int 	rowNum;					/* row number for master problem in solver */
}oneCut;

typedef struct {
	int    	cnt;                    /* number of cuts */
	oneCut  **vals;					/* values which define the set of cuts */
}cutsType;

typedef struct {
	int         k;                  /* number of iterations */
	int 		LPcnt; 				/* the number of LPs solved. */
    double		lb;					/* lower bound on cell objective function */
    int			lbType;				/* type of lower bound being used TRIVIAL if 0, else NONTRIVIAL */

    oneProblem  *master;            /* store master information */
	oneProblem 	*subprob;			/* store subproblem information */

	vector      candidX;            /* primal solution of the master problem */
	double      candidEst;          /* objective value master problem */

	vector      incumbX;			/* incumbent master solution */
	double      incumbEst;			/* estimate at incumbent solution */
	double 		quadScalar; 		/* the proximal parameter/quadratic scalar 'sigma' */
	BOOL        incumbChg;			/* set to be true if the incumbent solution has changed in an iteration */
	int         iCutIdx;			/* index of incumbent cut in cell->cuts structure */
	int         iCutUpdt;			/* iteration number when incumbent cut is updated */
	double      gamma;				/* improvement in objective function value */
	double      normDk_1;			/* (\Delta x^{k-1})^2 */
	double      normDk;				/* (\Delta x^k)^2 */

	vector      piS;                 /* subproblem dual information */
	double      mubBar;				/* dual slack information for subproblem */
	vector 		piM;				/* master dual information */
	vector      djM;                /* master reduced cost vector */

    int      	maxCuts;            /* maximum number of cuts to be used*/
	cutsType    *cuts;              /* optimality cuts */
	cutsType    *fcuts;             /* feasibility cuts */

	basisType	*basis;				/* hold unique basis identified */
	lambdaType 	*lambda;			/* holds dual solutions corresponding to rows effected by randomness */
	sigmaType 	*sigma;				/* holds $\pi \times \bar{b}$ and $\pi \times \bar{C} $ values */
	deltaType   *delta;				/* calculations based on realization and dual solutions observed */
	omegaType 	*omega;				/* all realizations observed during the algorithm */

    BOOL        optFlag;
	vector      pi_ratio;
    BOOL        dualStableFlag; 	/* indicates if dual variables are stable */

	int			feasCnt;			/* keeps track of the number of times infeasible candidate solution was encountered */
	BOOL		infeasIncumb;		/* indicates if the incumbent solution is infeasbible */
}cellType;

/* twoSD.c */
void parseCmdLine(string probName);
int readConfig(string inputDir);

/* algo.c */
int algo(oneProblem *orig, timeType *tim, stocType *stoc, string inputDir, string probName);
int solveCell(stocType *stoc, probType **prob, cellType *cell, string inputDir, string probName);
void writeStatistic(FILE **soln, probType **prob, cellType *cell, string probName, int numStages);
void cleanupAlgo(probType **prob, cellType *cell, int T);

/* setup.c */
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell);
cellType *newCell(stocType *stoc, probType **prob, vector xk);
void freeCellType(cellType *cell);

/* master.c */
int solveQPMaster(numType *num, sparseVector *dBar, cellType *cell, int IniRow, double lb);
int addCut2Master(cellType *cell, oneCut *cut, int lenX, double lb);
int constructQP(probType *prob, LPptr lp, vector incumbX);
int changeEtaCol(LPptr lp, int numRows, int numCols, int k, cutsType *cuts, double lb);
int updateRHS(LPptr lp, cutsType *cuts, int numIter, double lb);
int changeEtaCol(LPptr lp, int numRows, int numCols, int k, cutsType *cuts, double lb);
int updateRHS(LPptr lp, cutsType *cuts, int numIter, double lb);
int changeQPproximal(LPptr lp, int numCols, double sigma);
int changeQPrhs(probType *prob, cellType *cell);
int changeQPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector xk);
oneProblem *newMaster(oneProblem *orig, double lb);

/* cuts.c */
int formSDCut(probType *prob, cellType *cell, vector Xvect, int omegaIdx, BOOL newOmegaFlag, BOOL isIncumb);
oneCut *SDCut(numType *num, coordType *coord, basisType *basis, sigmaType *sigma, deltaType *delta, omegaType *omega, vector Xvect, int numSamples,
		BOOL *dualStableFlag, vector pi_ratio, double lb);
int computeIstar(numType *num, coordType *coord, basisType *basis, sigmaType *sigma, deltaType *delta, vector Xvect, vector PiCbarX, vector omegaVals, int obs,
		int numSamples, BOOL pi_eval, double *argmax, BOOL isNew);
oneCut *newCut(int numX, int numIstar, int numSamples);
cutsType *newCuts(int maxCuts);
int reduceCuts(cellType *cell, vector candidX, vector pi, int betaLen, double lb);
int dropCut(cellType *cell, int cutIdx);
double calcVari(double *x, double *mean_value, double *stdev_value, int batch_size);
void print_cut(cutsType *cuts, numType *num, int idx);
void freeOneCut(oneCut *cut);
void freeCutsType(cutsType *cuts);

/* subprob.c */
int solveSubprob(probType *prob, cellType *cell, vector Xvect, int omegaIdx, BOOL newOmegaFlag);
vector computeRHS(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, vector X, vector obs);
vector computeCostCoeff(numType *num, coordType *coord, sparseVector *dBar, vector obs, int offset);
void chgRHSwSoln(sparseVector *bBar, sparseMatrix *Cbar, vector rhs, vector X) ;
int chgRHSwObserv(LPptr lp, numType *num, coordType *coord, vector observ, vector spRHS, vector X);
int chgObjxwObserv(LPptr lp, vector cost, intvec indices, int rvdOmCnt, vector observ);
oneProblem *newSubproblem(oneProblem *subprob);

/* soln.c */
int checkImprovement(probType *prob, cellType *cell, int candidCut);
int replaceIncumbent(probType *prob, cellType *cell, double candidEst);
double maxCutHeight(cutsType *cuts, int currIter, vector xk, int betaLen, double lb);
double cutHeight(oneCut *cut, int currIter, vector xk, int betaLen, double lb);

/* optimal.c */
BOOL optimal(probType **prob, cellType *cell);
BOOL preTest(cellType *cell);
BOOL fullTest(probType **prob, cellType *cell);
cutsType *chooseCuts(cutsType *cuts, vector pi, int lenX);
void reformCuts(basisType *basis, sigmaType *sigma, deltaType *delta, omegaType *omega, numType *num, coordType *coord, cutsType *gCuts, int *observ, int k, int lbType, int lb, int lenX);
double calcBootstrpLB(probType *prob, vector incumbX, vector piM, vector djM, int currIter, double quadScalar, cutsType *cuts);
void empiricalDistribution(omegaType *omega, int *cdf);
void resampleOmega(intvec cdf, intvec observ, int numSamples);

/* evaluate.c */
int evaluate(FILE **soln, stocType *stoc, probType **prob, cellType *cell, vector Xvect);

/* stocUpdates.c */
int stochasticUpdates(cellType *cell, probType *prob, int omegaIdx, BOOL newOmegaFlag);

#endif /* TWOSD_H_ */
