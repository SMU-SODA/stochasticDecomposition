/*
 * algo.c
 *
 *  Created on: Apr 2, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "sdlp.h"

extern configType config;

int algo(oneProblem *orig, stocType *stoc, timeType *tim) {
	probType **prob = NULL;
	cellType **cell = NULL;
	vector	 observ;

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell) )
		goto TERMINATE;

	/* allocate memory to hold an observation */
	if ( !(observ = (vector) arr_alloc(stoc->numOmega, double)) )
		errMsg("allocation", "algo", "observ", 0);


	/* release memory allocated to different structures used in the algorithm */
	TERMINATE:
	cleanupAlgo(prob, cell, tim->numStages);
	mem_free(observ);

	return 0;
}//END algo()

void cleanupAlgo(probType **prob, cellType **cell, int T) {

	freeCellType(prob, cell, T);
	freeProbType(prob, T);

}//END cleanupAlgo()
