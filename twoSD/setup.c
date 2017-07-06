/*
 * setup.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"

extern configType config;

int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell) {
	vector	meanSol, lb;
	int 	t;

	/* setup mean value problem which will act as reference for all future computations */
	meanSol = meanProblem(orig, stoc);
	if ( meanSol == NULL ) {
		errMsg("setup", "setupAlgo", "failed to setup and solve mean value problem", 0);
		return 1;
	}

	/* calculate lower bounds for each stage */
	lb = calcLowerBound(orig, tim, stoc);
	if ( lb == NULL )  {
		errMsg("setup", "setupAlgo", "failed to compute lower bounds on stage problem", 0);
		return 1;
	}

	/* decompose the problem into master and subproblem */
	(*prob) = newProb(orig, stoc, tim, lb, config.TOLERANCE);
	if ( (*prob) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to update probType with elements specific to algorithm", 0);
		return 1;
	}

#ifdef DECOMPOSE_CHECK
	printDecomposeSummary(tim, (*prob));
#endif

	/* ensure that we have a linear programs at all stages */
	t = 0;
	while ( t < tim->numStages ) {
		if ( (*prob)[t++]->sp->type  != PROB_LP )
			printf("Warning :: Stage-%d problem is a mixed-integer program. Solving its linear relaxation.\n", t);
	}

	/* create the cells which will be used in the algorithms */
	(*cell) = newCell(stoc, (*prob), lb, meanSol);
	if ( (*cell) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to create the necessary cell structure", 0);
		return 1;
	}

	mem_free(meanSol); mem_free(lb);

	return 0;
}//END setupAlgo()

cellType *newCell(stocType *stoc, probType **prob, vector lb, vector meanSol) {
	cellType 	*cell;

	if ( !(cell = (cellType *) mem_malloc(sizeof(cellType))) )
		errMsg("allocation", "newCell", "cell", 0);
	cell->k = 0;
	cell->LP_cnt = 0;

	cell->cuts = new_cuts(prob->num->iter, p->num->mast_cols, 0);
	cell->lambda = new_lambda(length, 0, p->num->rv_rows, p->coord);
	cell->sigma = new_sigma(length, p->num->nz_cols, 0, p->coord);

	return cell;
}//END newCell

void freeCellType(probType *prob, cellType *cell) {

	if ( cell )
		mem_free(cell);

}//END freeCellType()
