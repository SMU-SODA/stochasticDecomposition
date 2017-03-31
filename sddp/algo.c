/*
 * algo.c
 *
 *  Created on: Mar 27, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *     Contact: harsha@smu.edu
 *
 */

#include "sddp.h"

int algo (oneProblem *orig, stocType *stoc, timeType *tim) {
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
