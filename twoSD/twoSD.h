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

typedef struct{
	double 	TOLERANCE; 			/* for zero identity test */
	int		MIN_ITER;			/* minimum number of iterations */
	int		MAX_ITER;			/* maximum number of iterations */
}configType;

typedef struct {
	double  alpha;                  /* scalar value for the righ-hand side */
	vector  beta;                   /* coefficients of the master problems's primal variables */
	int 	cutObs;					/* number of samples on which the given cut was based */
	int 	omegaCnt;				/* number of *distinct* observations on which the cut is based (this is also the length of istar) */
	intvec	iStar;					/* indices of maximal pi for each distint observation */
	BOOL	isIncumb;				/* indicates if the cut is an incumbent cut */
	double 	alphaIncumb;			/* right-hand side when using QP master, this is useful for quick updates */
	int 	slackCnt;				/* number of times a cut has been slack, used in deciding when the cut needs to be dropped */
	int 	rowNum;					/* row number for master problem in solver */
}oneCut;

typedef struct {
	int     cnt;                    /* number of cuts */
	oneCut  **val;
}cutsType;

/* To save time and space, Pi x b and Pi x C are calculated as soon as possible and stored in structures like sigma and delta.  Toward
 * this end, pixbCType represents a single calculation of pi X b (which is a scalar) and pi X C (which is a vector).*/
typedef struct{
	double 	b;
	vector 	C;
} pixbCType;

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
	intvec		ck; 				/* record the iteration # when sigma was created */
} sigmaType;

/* When calculating istar for a cut, it is useful to have two separate references into the sigma and delta structures, since each dual vector
 * is stored in two places -- part in sigma and part in delta.  The final entry in cut->istar[] will just be the _sigma_ field of this structure. */
typedef struct {
	int delta;
	int sigma;
} iType;

/* The delta matrix contains the values of lambda_pi X bOmega and lambda_pi X Comega for all values of pi and all observations of omega.
 * _col_ gives the column number of the each non-zero element in the multiplication of lambda_pi X Comega (the same elements are non-zero
 * each time).  _val_ is an array of vectors of (lambda_pi X bOmega, lambda_pi X Comega) pairs of calculations.  A row in _val_ corresponds
 * to a distinct dual vector, and a column in _val_ corresponds to a distinct observation of omega.  Thus, every pi-omega combination is
 * represented here, and the size of the delta matrix can be determine from lambda->cnt and omega->cnt.
 * Note that when elements of omega get dropped, vacant columns appear in delta.  This is ok, but be sure to loop carefully! */
typedef struct {
	pixbCType 	**vals;
} deltaType;

/**************************************************************************\
 ** Omega stores the set of observations which have been made so far.
 **
 **   Each observation consists of a vector of realizations of random
 ** variables with discrete distributions.  Since every distribution is
 ** discrete, an observation is just stored as a vector of indices into a
 ** distribution array containing the possible values.  _idx_ is an array
 ** of such vectors.  Each rv occurs in the R vector or T matrix which
 ** (along with the candidate X) make up the rhs of the subproblem.
 **
 **   The _row_ and _col_ arrays give the coordinates in R and T of each rv
 ** realization in a vector.  If _col_ is zero for a given entry, then the
 ** realization comes from R; otherwise, it comes from T.  The field _RT_
 ** represents a vector of actual realizations (as opposed to indices) of
 ** omega for both the Romega and Tomega structures (only one observation's
 ** worth of omega).
 **
 **   The _weight_ field specifies the number of times a particular outcome
 ** has been observed (it starts at 1 when the outcome is first generated,
 ** and increments every time the same outcome is observed again).  _cnt_
 ** just specifies the number of distinct outcomes which have been observed
 ** and stored in the omega structure.
 **
 **   If the problem gets too large, some unexciting omegas may be dropped.
 ** So, _filter_ (same length as _weight_) describes which vectors in the
 ** _idx_ array are actually filled with observations.  _next_ references
 ** the place to start in the _filter_ array when trying to find the next
 ** available position in _idx_.  _most_ represents the number of
 ** elements needed to store all elements in the _idx_ array, from
 ** the first to the last.  (It is the greatest index at which an
 ** observation is stored in the _idx_ array, plus one).   Finally,
 ** _last_ indicates the index of omega from which we started dropping
 ** omegas last time -- so we only "need" to go from omega.most down
 ** to omega.last when dropping them this time (we'll miss some).
 \**************************************************************************/
typedef struct {
	int 	cnt;
	intvec	weight;                 /* number of times that an omega is observed */
	vector	*vals;
} omegaType;

typedef struct{
    double		iterSolTime;        /* time of solving prob in each iter */
    double		totSolTime;         /* time of solving prob(cumulated) */
    double		iterCutGenTime;         /* time of generating cuts for each time */
    double  	totCutGenTime;      /* time of generating cuts(cumulated) */
    double		iterTime;        /* time for each iter (solveAgent) */
    double      totTime;     /* total time of doing a iteration for subprob */
//    vector  	masterIterTime;     /* time for each iter (within the while loop) */
//    double      totMasterIterTime;  /* time of cumulated iteration time */
}runTimeType;

typedef struct {
	int         k;                  /* number of iterations */
	oneProblem  *master;            /* store master information */
	oneProblem 	*subprob;			/* store subproblem information */
	vector      candidX;            /* primal solution of the master problem */
	double      candidEst;          /* objective value master problem */
    int         cCutIdx;            /* index of candidate cut in cell->cuts structure */
	vector      pi;                 /* dual slack information */
    vector      pi_ratio;
    BOOL        dualStableFlag; /* indicates if dual variables are stable */
	vector      di;                 /* reduced cost vector */
	double      mubBar;				/* dual slack information for subproblem */
	cutsType    *cuts;              /* optimality cuts */
	cutsType    *fcuts;             /* feasibility cuts */
	lambdaType 	*lambda;			/* holds dual solutions corresponding to rows effected by randomness */
	sigmaType 	*sigma;				/* holds $\pi \times \bar{b}$ and $\pi \times \bar{C} $ values */
	deltaType   *delta;				/* calculations based on realization and dual solutions observed */
	omegaType 	*omega;				/* all realizations observed during the algorithm */
	vector      incumbX;			/* incumbent master solution */
	double      incumbEst;			/* estimate at incumbent solution */
	double 		quadScalar; 		/* the proximal parameter/quadratic scalar 'sigma' */
	int 		LPcnt; 				/* the number of LPs solved. */
	int			feasCnt;			/* keeps track of the number of times infeasible candidate solution was encountered */
	BOOL		infeasIncumb;		/* indicates if the incumbent solution is infeasbible */
	BOOL		feasFlag;           /* indicates feasibility of a cell */
	double      incumbStdev;		/* standard deviation of incumbent estimate */
	BOOL        incumbChg;			/* set to be true if the incumbent solution has changed in an iteration */
	int         iCutIdx;			/* index of incumbent cut in cell->cuts structure */
	int         iCutUpdt;			/* iteration number when incumbent cut is updated */
	double      gamma;				/* improvement in objective function value */
	double      normDk_1;			/* (\Delta x^{k-1})^2 */
	double      normDk;				/* (\Delta x^k)^2 */
	double      optValM;            /* store the optimal value in the first stage */
    BOOL        optFlag;
	BOOL        newOmegaFlag;       /* set to true if a new omega is generated */
    int      	maxCuts;            /* maximum number of cuts to be used*/
    int			lbType;				/* type of lower bound being used TRIVIAL if 0, else NONTRIVIAL */
    vector      spRHS;              /* subproblem's rhs*/
    double      full_test_error;
    runTimeType *runTime;            /* record time */
    int         pushcnt;
}cellType;

/* twoSD.c */
void parseCmdLine(string probName);
int readConfig(string inputDir);

/* algo.c */
void printAlgoDetails(void *fptr);
void cleanupAlgo(probType **prob, cellType *cell, int T);

/* setup.c */
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell);
cellType *newCell(stocType *stoc, probType **prob, vector lb, vector meanSol);
void freeCellType(probType *prob, cellType *cell);

#endif /* TWOSD_H_ */
