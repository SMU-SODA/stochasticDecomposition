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
	int		cnt;

}cutsType;

typedef struct {

}lambdaType;

typedef struct {

}sigmaType;

typedef struct {
	int k;						/* iteration counter */
	oneProblem 	*master;		/* master problem structure */
	oneProblem 	*subprob;		/* subproblem structure */
	cutsType	*cuts;			/* minorants used to define lower bounding affine functions */
	lambdaType	*lambda;		/* holds relevant (for rows with random elements) duals */
	sigmaType	*sigma;			/* holds product information with respect to deterministic (mean/reference) problem */
	double 		quad_scalar; 	/* the quadratic scalar, 'sigma'. */
	int 		LP_cnt; 		/* number of LPs solved */
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
