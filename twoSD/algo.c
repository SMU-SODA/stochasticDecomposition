/*
 * algo.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"

extern string outputDir;
extern configType config;

int algo(string probName, oneProblem *orig, stocType *stoc, timeType *tim) {
	probType **prob = NULL;
	cellType *cell = NULL;
	vector	 observ;
	FILE *solnFile;

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell) )
		goto TERMINATE;

	/* allocate memory to hold an observation */
	if ( !(observ = (vector) arr_alloc(stoc->numOmega, double)) )
		errMsg("allocation", "algo", "observ", 0);

	printAlgoDetails(0);
	while (TRUE) {
#if VERBOSE
		printf("\nIteration %4d :: ", cell[0]->k+1); fflush(stdout);
#else
		if ( (cell->k) % 100 == 0)
			printf("\nIteration %4d :: ", cell->k); fflush(stdout);
#endif

//		/* if optimality conditions have been satisfied then break the while loop and exit. */
//		if ( optimal(prob, cell, tim->numStages) )
//			break;
//
//		/* generate a new observation */
//		generateOmega(stoc, observ, &config.RUN_SEED);
//
//		/* forward pass */
//		if (forwardPass(prob, cell, observ, tim->numStages)) {
//			errMsg("algorithm", "algo","failed in forward pass", 0);
//			goto TERMINATE;
//		}
//
//		/* backward pass */
//		if (backwardPass(prob, cell, observ, tim->numStages)) {
//			errMsg("algorithm", "algo", "failed in backward pass", 0);
//			goto TERMINATE;
//		}
//
//		/* check to see if there is an improvement */
//		checkImprovement(prob, cell, tim->numStages);
	}

	printf("\n\n-------------------------------------------------------------------------------------------------------------------------------\n");
	printf("Successfully completed excecution of SDLP algorithm on %s.\n", probName);
	printf("-------------------------------------------------------------------------------------------------------------------------------\n");

	/* Print solution details */
	solnFile = openFile(outputDir, "detailedSDLPsols.dat", "w");
	/* release memory allocated to different structures used in the algorithm */
	TERMINATE:
	cleanupAlgo(prob, cell, tim->numStages);
	mem_free(observ);

	return 0;
}//END algo()

void printAlgoDetails(void *fptr) {

	fprintf(fptr, "\n-------------------------------------------------------------------------------------------------------------------------------\n");
	fprintf(fptr, "Starting two-stage stochastic decomposition algorithm\n");
	fprintf(fptr, "Minimum iterations: %d\n", config.MIN_ITER);
	fprintf(fptr, "-------------------------------------------------------------------------------------------------------------------------------\n");

}//END printAlgoDetails()

void cleanupAlgo(probType **prob, cellType *cell, int T) {

	freeCellType(prob, cell);
	freeProbType(prob, T);

}//END cleanupAlgo()
