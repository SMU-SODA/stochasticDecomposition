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

#endif /* SDLP_H_ */
