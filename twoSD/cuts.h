/*
 * cuts.h
 *
 *  Created on: Sep 11, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#ifndef CUTS_H_
#define CUTS_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"

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

oneCut *newCut(int numX, int numIstar, int numSamples);
cutsType *newCuts(int maxCuts);
double maxCutHeight(cutsType *cuts, int currIter, vector xk, int betaLen, double lb);
double cutHeight(oneCut *cut, int currIter, vector xk, int betaLen, double lb);
int reduceCuts(oneProblem *master, cutsType *cuts, int *iCutIdx, vector candidX, vector pi, int betaLen, double lb, int currentIter, double TOLERANCE);
int dropCut(oneProblem *master, cutsType *cuts, int cutIdx, int *iCutIdx, int currentIter);
void freeOneCut(oneCut *cut);
void freeCutsType(cutsType *cuts);


#endif /* CUTS_H_ */
