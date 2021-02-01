/*
 * twoSD.h
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
  *  
 *
 *   Edited on: Dec 1, 2020 as part of the SD-integer project
 *      Author: Siavash Tabrizian
 * Institution: Southern Methodist University
 *
 * Please send you comments or bug report to stabrizian (at) smu (dot) edu
 *
 */

#ifndef TWOSD_H_
#define TWOSD_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"
#include "stoc.h"

#undef experiment

#define TRIVIAL 0
#define NONTRIVIAL 1

#undef ALGO_CHECK
#undef BNC_CHECK
#undef BATCH_CHECK
#define clean_master /* clean the master problem at the begining of each node */

#undef LPMIP_PRINT
#undef PHASE1ANLYS

#undef CALLBACK_CHECK
#define CALLBACK_WRITE_LP

#undef UserMIPcutsActive

#if defined(UserMIPcutsActive)
#define GMIcutsActive
#define GMIstreghtcuts1Active
#define GMIstreghtcuts2Active
#define MIRcutsActive
#define MIRSubbaddActive
#else
#undef CpxGMICutsActive
#undef CpxMIRCutsActive
#endif // defined(UserMIPcutsActive)

enum smipSolver {
	MILP,
	CUSTOM,
	CALLBACK
};

enum cutForm {
	Regular,
	Callback,
	Feasible,
	MIR,
	GMI
};

/* A data structure which holds on the configuration information about the algorithm. Most of these configuration parameters are read from a
-configuration file. These elements, once set during initialization, are not modified during the course of the algorithm. */
typedef struct{
	int		NUM_REPS;			/* Maximum number of replications that can be carried out. */
	long long *RUN_SEED;		/* seed used during optimization */
	double 	TOLERANCE; 			/* for zero identity test */
	int		MIN_ITER;			/* minimum number of iterations */
	int		MAX_ITER;			/* maximum number of iterations */
	int     MAX_NODES;          /* Maximum number of visited nodes */
	int     VAR_STR;            /* variable selection strategy */
	int     BRN_STR;            /* branching strategy */
	int		MAX_ITER_CLBK;		/* maximum number of iterations inside the callback */
	int		MASTER_TYPE;		/* type of master problem */
	int		TAU;				/* Frequency at which the incumbent is updated */
	double	MIN_QUAD_SCALAR;	/* Minimum value for regularizing parameter */
	double	EPSILON;			/* Optimality gap */

	int		EVAL_FLAG;
	int		Pi_EVAL_FLAG;       /* If this flag is on only the recent duals will be evaluated in the B&B tree */
	int     HEURST_FLAG;        /* using heuristic in bnc */
	int		NUM_EVALS;
	long long *EVAL_SEED;
	int		EVAL_MIN_ITER;
	double	EVAL_ERROR;

	int		CUT_MULT;			/* Determines the number of cuts to be used for approximate */
	double 	MAX_QUAD_SCALAR;	/* Maximum value for regularizing parameter */
	double 	MIN_X;	            /* Mimimum X norm difference for incumbent update */
	double	R1;
	double	R2;
	double	R3;
	int		DUAL_STABILITY;		/* Determine whether dual stability is to be checked. */
	int		PI_EVAL_START;		/* The minimum number of samples before dual stability test is conducted. */
	int		PI_CYCLE;			/* Frequency of updating the dual stability ratio */
	int		BOOTSTRAP_REP;		/* Number of boot-strap replications in full optimality test */
	double	PERCENT_PASS;		/* percentage of bootstrap replications need to be satisfied */
	int		SCAN_LEN;			/* window size over which the stability of dual vertex set is measured.*/
	double  PRE_EPSILON;		/* gap used for preliminary optimality test */

	int 	SMIP;
	int     NodeNum;
	int 	ALGO;
	double 	SMIP_OPTGAP;
	int 	MULTIPLE_REP;		/* When multiple replications are needed, set this to (M), else (0) */
	int		COMPROMISE_PROB;	/* Compromise solution created and solved for compromise solution. */

	int 	SAMPLE_INCREMENT;	/* Number of new observations added to the sample */
	double     INT_TOLERANCE;      /* Tolerence for call a variable integer*/
}configType;

typedef struct {
	bool    isAvctive;				/* is the cut active in the master problem */
	double  alpha;                  /* scalar value for the righ-hand side */
	dVector  beta;                  /* coefficients of the master problems's primal variables */
	int 	numSamples;				/* number of samples on which the given cut was based */
	int 	omegaCnt;				/* number of *distinct* observations on which the cut is based (this is also the length of istar) */
	iVector	iStar;					/* indices of maximal pi for each distint observation */
	bool	isIncumb;				/* indicates if the cut is an incumbent cut */
	double 	alphaIncumb;			/* right-hand side when using QP master, this is useful for quick updates */
	int 	slackCnt;				/* number of times a cut has been slack, used in deciding when the cut needs to be dropped */
	int 	rowNum;					/* row number for master problem in solver */
	cString	name;
}oneCut;

typedef struct {
	int     cnt;                    /* number of cuts */
	oneCut  **vals;					/* A vector of oneCut type, each element has information about a particular cut. */
}cutsType;

typedef struct {
	double	repTime;
	double 	iterTime;
	double 	masterIter;
	double 	subprobIter;
	double 	optTestIter;
	double 	argmaxIter;
	double 	iterAccumTime;
	double 	masterAccumTime;
	double 	subprobAccumTime;
	double 	optTestAccumTime;
	double 	argmaxAccumTime;
}runTime;

typedef struct {
	int seqid;                      /* sequence id of the node in the B&B tree work with CPXgetcallbackseqinfo*/
	int nodeidx;                    /* index of the node when the callback is called with CPXgetcallbacknodeinfo*/
	int nodeNum;
	int sol_size;
	dVector sol;					/* candidate solution of node */
	double LB;
	double UB;
	int    omegaSize;               /* Store the sample size associate with this node */
	bool   isInt;					/* is the candid solution integer */
	dVector 	piM;				/* master dual information */
	dVector     djM;                /* master reduced cost dVector */
	int    mar;
}nodeInfo;

// info that can be passed from the bnc
typedef struct {
	iVector disjnct;				/* shows which variables are disjuncted */
	dVector vals;                   /* value of variables including those are disjuncted */
}bncInfoSummary;

typedef struct {
	cString probName;                /* problem name */
									 /* EXperiment results */
	int         ssize;
	double      ssize_var;
	double      ssize_stdev;
	double      LBest;
	double      LBstdev;
	double      LBvar;
	double      UBest;
	double      UBstdev;
	double      UBvar;
	double      tot_time;
	double      tot_time_var;
	double      tot_time_stdev;
}expSummary;

typedef struct {
	int         k;                  /* number of iterations */
	int         ki;                 /* number of iterations cumulatively after callback*/
	int         kii;                /* number of iterations after callback*/
	int         gk;                 /* number of GMI cuts */
	int         mk;                 /* number of MIR cuts */
	int 		sampleSize;			/* total number of observations currently being used, that is the sample size. */
	int 		LPcnt; 				/* the number of LPs solved. */
    double		lb;					/* lower bound on cell objective function */
    int			lbType;				/* type of lower bound being used TRIVIAL if 0, else NONTRIVIAL */
	bool		callback;			/* flag to indicate if the cell is being solved in callback phase (true) */
	double      meanVal;            /* optimal objective of the mean value problem */

    oneProblem  *master;            /* store master information */
	oneProblem 	*subprob;			/* store subproblem information */

	oneProblem 	*PHsubprob;			/* store Progressive Hedging subproblem information */

	dVector      candidX;           /* primal solution of the master problem */
	double      candidEst;          /* objective value master problem */

	dVector     incumbX;			/* incumbent master solution */
	dVector     incumbMIPX;			/* incumbent master solution after MIP procedure */
	double      incumbEst;			/* estimate at incumbent solution */
	double 		quadScalar; 		/* the proximal parameter/quadratic scalar 'sigma' */
	bool        incumbChg;			/* set to be true if the incumbent solution has changed in an iteration */
	int         iCutIdx;			/* index of incumbent cut in cell->cuts structure */
	int         iGCutIdx;			/* index of GMI incumbent cut in cell->MIPcuts structure */
	int         iMCutIdx;			/* index of MIR incumbent cut in cell->MIPcuts structure */
	int         iCutUpdt;			/* iteration number when incumbent cut is updated */
	double      gamma;				/* improvement in objective function value */
	double      normDk_1;			/* (\Delta x^{k-1})^2 */
	double      normDk;				/* (\Delta x^k)^2 */

	dVector 	piM;				/* master dual information */
	dVector     djM;                /* master reduced cost dVector */

    int      	maxCuts;            /* maximum number of cuts to be used*/
	int      	maxMIPCuts;         /* maximum number of MIP cuts to be used*/
	int         etaIdx;             /* Index of eta column after starting the B&B*/

	cutsType    *cuts;              /* optimality cuts */
	cutsType    *MIRcuts;           /* MIP feasibility cuts */
	cutsType    *GMIcuts;           /* MIP feasibility cuts */
	cutsType    *fcuts;             /* feasibility cuts */
    cutsType 	*fcutsPool;			/* Pool of feasibility cuts */
    int			fUpdt[2];			/* coordinate in delta structure for which the updates have been carried out */

	lambdaType 	*lambda;			/* holds dual solutions corresponding to rows effected by randomness */
	sigmaType 	*sigma;				/* holds $\pi \times \bar{b}$ and $\pi \times \bar{C} $ values */
	deltaType   *delta;				/* calculations based on realization and dual solutions observed */
	omegaType 	*omega;				/* all realizations observed during the algorithm */
	basisType	*basis;				/* hold unique basis identified */

    bool        optFlag;			/* Optimality flag */
	bool        MIPFlag;			/* MIP feasibility flag */
	dVector     pi_ratio;			/* Pi ratios over a window of selected size (determined by tolerance level) */
    bool        dualStableFlag; 	/* indicates if dual variables are stable */

    bool 		optMode;			/* When false, the algorithm tries to resolve infeasibility */

    bool		spFeasFlag;			/* Indicates whether the subproblem is feasible */
	int			feasCnt;			/* keeps track of the number of times infeasible candidate solution was encountered */
	bool		infeasIncumb;		/* indicates if the incumbent solution is infeasible */

	runTime		time;				/* Run time structure */

	///// B&B paramters
	bncInfoSummary *bncInfo;
	bool isinBnB;
	int rownum;
	sampleType	*sample;
	nodeInfo    *nodeSol;			/* solution of nodes that are discovered in B&B */
	char        **cur_rowname;		/* row names in BnB */

}cellType;

typedef struct {
	oneProblem	*sp;				/* compromise problem */
	int 		cnt;				/* number of replications */
	iVector 	ck;					/* number of iterations for each replication */
	dVector		objLB;				/* replication lower bound */
	dVector		objUB;				/* replication upper bound, if batch solution is evaluated */
	double		objComp;			/* optimal value of compromise problem */
	double		quadScalar;			/* average proximal terms */
	dVector		*incumbX;			/* batch incumbent solution */
	dVector		compromiseX;		/* compromise solution */
	dVector		avgX;				/* average solution across batches */
}batchSummary;


/// BnC structures 
/* A data structure which holds on the node informations in the branch and bound tree */
struct BnCnodeType {
	int		key;			             /* key id of the node. */
	int		parentkey;			         /* parent key id of the node. */
	int		depth;			             /* depth of the node on the tree. */
	int     varId;                       /* index of the variable that is disjuncted in this node */
	int     stInt;                       /* index of the first integer variable */
	int     edInt;                       /* index of the last integer variable */
	int     parentnumSamp;               /* sample size of the parent */
	int     numSamp;                     /* sample size of the current node */
	double  fracPi;                      /* (\tilde(Pi)_{n-1} + \Delta\Pi_{n})/\Pi_{n} */
	int     tightPi;                     /* number of pis used in tight cuts */
	int     Lambdasize;                  /* total lamda size */
	int     partightPi;                  /* number of pis used in tight cuts of th parent node */
	int     parLambdasize;               /* total lamda size the parent node */
	int     parparinit;                  /* initial lambda loop for the next iterations from parent of the parent */
	double  fracVal;                     /* disjuncted fractional value for the varId */
	double  parobjVal;                   /* objective value of the parent node */
	iVector  disjncs;                    /* list of distjuctive cuts on variables 0: not added 1: disjnct is added */
	dVector  * disjncsVal;               /* list of distjuctive cuts values on variables - it has upper and lower limits */
	iVector	IncumbiStar;				 /* indices of maximal pi for each distint observation for incumbent cuts */
	iVector	ParIncumbiStar;				 /* indices of maximal pi for each distint observation for incumbent cuts of the parent */
	int numVar;                          /* number of variables */
	int numRows;
	double LB;                           /* lower bounds and upper bounds at this node */
	double UB;                           /* lower bounds and upper bounds at this node */
	double objval;                       /* objective value of the node problem  */
	dVector vars;                        /* the values of variables */
	dVector duals;                       /* the values of the duals of the node problem */
	bool   isInt;                        /* is the solution obtained from the node integer */
	bool   isActive;                     /* the node is active or not */
	bool   isSPopt;                      /* the estimation of the node completed */
	bool   isfathomed;                   /* the node is fathomed or not */
	bool   isleft;                       /* the node is in the left of its parent */
	bool   ishrstic;                     /* the node is a heuristic node */
	struct BnCnodeType *left, *right;    /* pointer to the left and right children. */
	struct BnCnodeType *nextnode;        /* next node in the queue. */
	struct BnCnodeType *prevnode;        /* next node in the queue. */
};


////// Subroutines
/* algo.c */
int algo(oneProblem *orig, timeType *tim, stocType *stoc, cString inputDir, cString probName);
int solveCell(stocType *stoc, probType **prob, cellType *cell);
int QPtoLP(stocType *stoc, probType **prob, cellType *cell, int toMIP);
int LPtoMILP(stocType *stoc, probType **prob, cellType *cell);
int mainloopSDCell(stocType *stoc, probType **prob, cellType *cell, bool *breakLoop, dVector observ);
void writeOptimizationSummary(FILE *soln, FILE *incumb, probType **prob, cellType *cell, bool header);
int phase_one_analysis(stocType *stoc, probType **prob, cellType *cell);
void printNodeInfo(nodeInfo    *nodeSol, int Nodecnt);
void getRowNameMaster(cellType *cell);

/* setup.c */
int readConfig(cString path2config, cString inputDir);
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell, batchSummary **batch, dVector *meanSol, dVector lb);
int setupClone(oneProblem *orig, stocType *stoc, timeType *tim, probType ***cloneprob, cellType **clonecell, dVector *meanSol, dVector lb);
cellType *newCell(stocType *stoc, probType **prob, dVector xk);
int cleanCellType(cellType *cell, probType *prob, dVector xk);
int cleanBnCCellType(cellType *cell, probType *prob, dVector xk);
void freeCellType(cellType *cell);

/* master.c */
int solveQPMaster(numType *num, sparseVector *dBar, cellType *cell, double lb);
int solveLPMaster(numType *num, sparseVector *dBar, cellType *cell, double lb);
int addCut2Master(oneProblem *master, oneCut *cut, dVector vectX, int lenX);
int addMIPCut2Master(oneProblem *master, oneCut *cut, dVector vectX, int lenX, bool GMI);
int constructQP(probType *prob, cellType *cell, dVector incumbX, double quadScalar);
int changeEtaColMIP(LPptr lp, int numRows, int numCols, int currSampleSize, cutsType *SDcuts, cutsType *MIRcuts, cutsType *GMIcuts, int iter, char **cur_rowname);
int changeEtaCol(LPptr lp, int numRows, int numCols, int currSampleSize, cutsType *cuts);
int updateRHS(LPptr lp, cutsType *cuts, int numIter, double lb);
int changeQPproximal(LPptr lp, int numCols, double sigma);
int changeQPrhs(probType *prob, cellType *cell, dVector xk);
int revchangeQPrhs(probType *prob, cellType *cell, dVector xk);
int changeQPbds(LPptr lp, int numCols, dVector bdl, dVector bdu, dVector xk, int offset);
oneProblem *newMaster(oneProblem *orig, double lb);
int recoverX(dVector sol, dVector incumb, dVector candid, iVector disjnct, dVector vals, int cnt);

/* cuts.c */
int formSDCut(probType **prob, cellType *cell, dVector Xvect, double lb, int inCallback);
int formGMICut(probType **prob, cellType *cell, dVector Xvect, double lb);
int formMIRCut(probType **prob, cellType *cell, dVector Xvect, double lb);
oneCut **pureMIRCut(probType **prob, cellType *cell, dVector Xvect, double lb);
oneCut **purefracGMICut(probType **prob, cellType *cell, dVector Xvect, double lb);
oneCut *SDCut(numType *num, coordType *coord, basisType *basis, sigmaType *sigma, deltaType *delta, omegaType *omega, sampleType *sample,
		dVector Xvect, int numSamples, bool *dualStableFlag, dVector pi_ratio, int numIter, double lb);
oneCut *newCut(int numX, int numIstar, int numSamples);
cutsType *newCuts(int maxCuts);
int reduceCuts(cellType *cell, dVector candidX, dVector pi, int betaLen, double lb);
int dropCut(cellType *cell, int cutIdx);
double calcVariance(double *x, double *mean_value, double *stdev_value, int batch_size);
void printCut(oneCut *cut, int betaLen);
void freeOneCut(oneCut *cut);
void freeCutsType(cutsType *cuts, bool partial);
double calc_var(double *x, double *mean_value, double *stdev_value, int batch_size);

/* soln.c */
int checkImprovement(probType *prob, cellType *cell, int candidCut);
int checkImprovement_callback(probType *prob, cellType *cell, int candidCut);
int replaceIncumbent(probType *prob, cellType *cell, double candidEst);
double maxCutHeight(cutsType *cuts, int currIter, dVector xk, int betaLen, double lb);
double cutHeight(oneCut *cut, int currIter, dVector xk, int betaLen, double lb);

/* optimal.c */
bool optimal(probType **prob, cellType *cell);
bool LPoptimal(probType **prob, cellType *cell);
bool preTest(cellType *cell);
bool fullTest(probType **prob, cellType *cell);
cutsType *chooseCuts(cutsType *cuts, dVector pi, int lenX);
void reformCuts(basisType *basis, sigmaType *sigma, deltaType *delta, omegaType *omega, numType *num, coordType *coord,
		cutsType *gCuts, int *observ, int sampleSize, int lbType, int lb, int lenX);
double calcBootstrpLB(probType *prob, dVector incumbX, dVector piM, dVector djM, int currIter, double quadScalar, cutsType *cuts, dVector bl, dVector bu);
void empiricalDistribution(omegaType *omega, int *cdf);
void resampleOmega(iVector cdf, iVector observ, int numSamples);
bool IPoptimal(probType **prob, cellType *cell);
bool IPpreTest(cellType *cell);

/* compromise.c */
int buildCompromise(probType *prob, cellType *cell, batchSummary *batch);
int solveCompromise(probType *prob, batchSummary *batch);
int addBatchEquality (probType *prob, batchSummary *batch);
batchSummary *newBatchSummary(probType *prob, int numBatch);
void freeBatchType(batchSummary *batch);

/* evaluate.c */
int evaluate(FILE *solnSumm,FILE *soln, stocType *stoc, probType **prob, oneProblem *subprob, dVector Xvect);
void writeEvaluationSummary(FILE *soln, double mean, double stdev, int cnt);

/* inout.c */
void writeOptimizationStatistics(FILE *soln, FILE *incumb, probType **prob, cellType *cell, int rep);
void writeEvaluationStatistics(FILE *soln, double mean, double stdev, int cnt);
void writeOmegaHist(FILE *soln, cellType *cell, int rep);
void printOptimizationSummary(cellType *cell);
void printEvaluationSummary(FILE *soln, double mean, double stdev, int cnt);

/* bnc.c */
int sumintVec(iVector a, int len);
bool isInVec(iVector vec, int len, int val);
struct BnCnodeType *newrootNode(int numVar, double LB, double UB, oneProblem * orig);
struct BnCnodeType *newNode(int key, struct BnCnodeType * parent, double fracVal, int varId, bool isleft);
int addBnCDisjnct(cellType *cell, dVector  *disjncsVal, int numCols, struct BnCnodeType * node);
double solveNode(stocType *stoc, probType **prob, cellType *cell, struct BnCnodeType *node, cString pname);
int          freeNodes(struct BnCnodeType *root);
int 	     freeNode(struct BnCnodeType *node);
struct BnCnodeType *copyNode(struct BnCnodeType *node, double thresh);
struct BnCnodeType *nextNode(struct BnCnodeType *node);
struct BnCnodeType     *successorNode(struct BnCnodeType *node);
struct BnCnodeType *insertNode(struct BnCnodeType *node, struct BnCnodeType *activenode);
int branchbound(stocType *stoc, probType **prob, cellType *cell, double LB, double UB);
int branchVar(struct BnCnodeType *node, int strategy);
int branchNode(stocType *stoc, probType **prob, cellType *cell, struct BnCnodeType *node, struct BnCnodeType **activeNode);
int getfirstLeaf(int depth);
int getnodeIdx(int depth, int key, int isleft);
void fracLamda(cellType *cell, struct BnCnodeType *node);
void truncate(dVector var, dVector lb, dVector ub, int cnt);

#endif /* TWOSD_H_ */
