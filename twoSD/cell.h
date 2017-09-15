/*
 * cell.h
 *
 *  Created on: Sep 12, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#ifndef CELL_H_
#define CELL_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"
#include "stoc.h"

#define TRIVIAL 0
#define NONTRIVIAL 1

/* A data structure which holds on the configuration information about the algorithm. Most of these configuration parameters are read from a
configuration file. These elements, once set during initialization, are not modified during the course of the algorithm. */
typedef struct{
	long long RUN_SEED;			/* seed used during optimization */
	double 	TOLERANCE; 			/* for zero identity test */
	int		MIN_ITER;			/* minimum number of iterations */
	int		MAX_ITER;			/* maximum number of iterations */
	int		MASTER_TYPE;		/* type of master problem */
	int		TAU;				/* Frequency at which the incumbent is updated */
	double	MIN_QUAD_SCALAR;	/* Minimum value for regularizing parameter */
	double	EPSILON;			/* Optimality gap */

	int		EVAL_FLAG;
	long long EVAL_SEED;
	int		EVAL_MIN_ITER;
	double	EVAL_ERROR;

	int		CUT_MULT;			/* Determines the number of cuts to be used for approximate */
	double 	MAX_QUAD_SCALAR;	/* Maximum value for regularizing parameter */
	double	R1;
	double	R2;
	double	R3;
	int		PI_EVAL_START;
	int		PI_CYCLE;
	int		BOOTSTRAP_REP;		/* Number of boot-strap replications in full optimality test */
	double	PERCENT_PASS;		/* percentage of bootstrap replications need to be satisfied */
	int		SCAN_LEN;			/* window size over which the stability of dual vertex set is measured.*/
	double  PRE_EPSILON;		/* gap used for preliminary optimality test */

	int		MAX_OBS;			/* Maximum number of iterations before which SAA is invoked */
	double  SUBPROB_SAMPLE_PCT;		/* Fraction of subproblem being solved in an iteration */
	long long SUBPROB_SAMPLE_SEED;	/* Seed used to sample the subproblems. */
	int		SAA; 				/* Use SAA when continuous distribution in stoch file (1), or not (0) */
}configType;

/* The oneCut and cutsType data structures will be used to hold all information which can completely define the affine minorants (cuts) which
 * are used to compute the lower bounding function approximations */
typedef struct {
	double  alpha;                  /* scalar value for the right-hand side */
	vector  beta;                   /* coefficients of the master problems's primal variables */
	int 	numSamples;				/* number of samples on which the given cut was based */
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
	double		ub;					/* upper bound */

	vector      incumbX;			/* incumbent master solution */
	double      incumbEst;			/* estimate at incumbent solution */
	double 		quadScalar; 		/* the proximal parameter/quadratic scalar 'sigma' */
	BOOL        incumbChg;			/* set to be true if the incumbent solution has changed in an iteration */
	int         iCutIdx;			/* index of incumbent cut in cell->cuts structure */
	int         iCutUpdt;			/* iteration number when incumbent cut is updated */
	double      gamma;				/* improvement in objective function value */
	double      normDk_1;			/* (\Delta x^{k-1})^2 */
	double      normDk;				/* (\Delta x^k)^2 */

	vector 		piM;				/* master dual information */
	vector      djM;                /* master reduced cost vector */

    int      	maxCuts;            /* maximum number of cuts to be used*/
	cutsType    *cuts;              /* optimality cuts */
	cutsType    *fCuts;             /* feasibility cuts */
    cutsType 	*fCutsPool;			/* Pool of feasibility cuts */
    int			fUpdt[2];			/* coordinate in delta structure for which the updates have been carried out */

	basisType	*basis;				/* hold unique basis identified */
	lambdaType 	*lambda;			/* holds dual solutions corresponding to rows effected by randomness */
	sigmaType 	*sigma;				/* holds $\pi \times \bar{b}$ and $\pi \times \bar{C} $ values */
	deltaType   *delta;				/* calculations based on realization and dual solutions observed */
	omegaType 	*omega;				/* all realizations observed during the algorithm */

    BOOL        optFlag;
	vector      pi_ratio;
    BOOL        dualStableFlag; 	/* indicates if dual variables are stable */

    BOOL		optMode;
    BOOL		spFeasFlag;			/* Indicates whether the subproblem is feasible */
    int			feasCnt;			/* keeps track of the number of times infeasible candidate solution was encountered */
	BOOL		infeasIncumb;		/* indicates if the incumbent solution is infeasbible */
}cellType;

/* setup.c */
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell);
cellType *newCell(stocType *stoc, probType **prob, vector xk);
oneProblem *newMaster(oneProblem *orig, double lb);
int constructQP(probType *prob, cellType *cell, vector incumbX, double quadScalar);
int changeQPproximal(LPptr lp, int numCols, double sigma);
int changeQPrhs(probType *prob, cellType *cell, vector xk);
int changeQPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector xk);
int readConfig();
void freeCellType(cellType *cell);

/* cuts.c */
int addCut2Master(cellType *cell, oneCut *cut, BOOL scaleCut, int lenX, double lb);
int replaceIncumbent(probType *prob, cellType *cell, double candidEst);;
oneCut *newCut(int numX, int numIstar, int numSamples);
cutsType *newCuts(int maxCuts);
double maxCutHeight(cutsType *cuts, vector xk, int betaLen, BOOL scaleCut, int currIter, double lb);
double cutHeight(oneCut *cut, vector xk, int betaLen, BOOL scaleCut, int currIter, double lb);
int reduceCuts(oneProblem *master, cutsType *cuts, BOOL scaleCut, vector candidX, vector pi, int betaLen, double lb, int currentIter, int *iCutIdx, double TOLERANCE);
int dropCut(oneProblem *master, cutsType *cuts, int cutIdx, int *iCutIdx);
void freeOneCut(oneCut *cut);
void freeCutsType(cutsType *cuts);

int formFeasCut(probType *prob, cellType *cell, BOOL *newOmegaFlag, BOOL newBasisFlag);
int updtFeasCutPool(numType *num, coordType *coord, cutsType *fCutsPool, int fUpdt[2], basisType *basis, sigmaType *sigma, deltaType *delta, omegaType *omega,
		BOOL newOmegaFlag, BOOL newBasisFlag, int currentIter);
int add2CutPool(cutsType *cuts, double alpha, vector beta, int betaLen, int numOmega, int numSamples);
int checkFeasCutPool(cutsType *cutPool, cutsType *cutsAdded, int betaLen, vector incumbX, vector candidX, BOOL *infeasIncumb);
int addfCut2Master(LPptr lp, oneCut *cut, vector incumbX, int lenX, int optCuts, int idx);

#endif /* CELL_H_ */
