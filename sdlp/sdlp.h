/*
 * sdlp.h
 *
 *  Created on: Apr 2, 2017
 *      Author: gjharsha
 */

#ifndef SDLP_H_
#define SDLP_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"

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
	int			k;

}cellType;

/* sdlp.c */
void parseCmdLine(string probName);
int readConfig(string inputDir);

/* algo.c */
void cleanupAlgo(probType **prob, cellType **cell, int T);

/* setup.c */
int algo(oneProblem *orig, stocType *stoc, timeType *tim);
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType ***cell);

#endif /* SDLP_H_ */
