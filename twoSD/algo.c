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

int algo(oneProblem *orig, timeType *tim, stocType *stoc, string inputDir, string probName) {
	vector	 xk = NULL, lb = NULL;
	probType **prob = NULL;
	cellType *cell = NULL;
	double   totRunTime;

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell) )
		goto TERMINATE;

	/* Use two-stage algorithm to solve the problem */
	if ( solveCell(stoc, prob, cell, inputDir, probName) ) {
		errMsg("algorithm", "algo", "failed to solve the cells using MASP algorithm", 0);
		goto TERMINATE;
	}

	/* Write the solutions statistics */
	writeStat(prob[0], cell, totRunTime, probName);

	/* free up memory before leaving */
	if (xk) mem_free(xk);
	if (lb) mem_free(lb);
	freeCellType(cell);
	freeProbType(prob, 2);
	return 0;

	TERMINATE:
	if(xk) mem_free(xk);
	if(lb) mem_free(lb);
	if(cell) freeCellType(cell);
	if(prob) freeProbType(prob, 2);
	return 1;
}//END algo()

int solveCell(stocType *stoc, probType **prob, cellType *cell, string inputDir, string probName) {
	vector 	observ;
	int		m, omegaIdx;
	BOOL 	newOmegaFlag;

	/* -+-+-+-+-+-+-+-+-+-+-+-+-+-+- Main Algorithm -+-+-+-+-+-+-+-+-+-+-+-+-+-+- */
	if ( !(observ = (vector) arr_alloc(stoc->numOmega + 1, double)) )
		errMsg("allocation", "solveMASP", "observ", 0);

	/******* 0. Initialization: The algorithm begins by solving the master problem as a QP *******/
	while (cell->optFlag == FALSE && cell->k < config.MAX_ITER) {
		cell->k++;

#if 0
		printf("\nIteration-%d :: \n", cell->k);
#else
		if ( (cell->k -1) % 50 == 0)
			printf("\nIteration-%4d: ", cell->k);
#endif

		/******* 1. Optimality tests *******/
		if (optimal(prob, cell))
			break;

		/******* 2. Generate new observation, and add it to the set of observations *******/
		/* (a) Use the stoc file to generate observations */
		generateOmega(stoc, observ, &config.RUN_SEED);

		/* (b) Since the problem already has the mean values on the right-hand side, remove it from the original observation */
		for ( m = 0; m < stoc->numOmega; m++ )
			observ[m] -= stoc->mean[m];

		/* (d) update omegaType with the latest observation. If solving with incumbent then this update has already been processed. */
		omegaIdx = calcOmega(observ - 1, 0, prob[1]->num->numRV, cell->omega, &newOmegaFlag);

		/******* 3. Solve the subproblem with candidate solution, form and update the candidate cut *******/
		if ( formSDCut(prob[1], cell, cell->candidX, observ, TRUE) ) {
			errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
			return 1;
		}
		/******* 4. Solve subproblem with incumbent solution, and form an incumbent cut *******/
		if (((cell->k - cell->iCutUpdt) % config.TAU == 0 ) ) {
			if ( formSDCut(prob[1], cell, cell->incumbX, observ, FALSE) ) {
				errMsg("algorithm", "solveCell", "failed to create the incumbent cut", 0);
				return 1;
			}
		}

		/******* 5. Check improvement in predicted values at candidate solution *******/
		if ( !(cell->incumbChg) && cell->k > 1) {
			/* If the incumbent has not changed in the current iteration */
			checkImprovement(prob[0], cell);

			if (cell->incumbChg)
				/******* 6. If the incumbent solution has changed then update all things concerned with the incumbent */
				if ( constructQP(cell->master->lp, prob[0]->num->cols, cell->quadScalar) ) {
					errMsg("algorithm", "solveMASP", "failed to change the proximal term", 0);
					return 1;
				}
		}

		/******* 6. Solve the master problem to obtain the new candidate solution */
		if ( solveQPMaster(prob[0]->num, prob[0]->dBar, cell, prob[0]->sp->mar, prob[0]->lb) ) {
			errMsg("algorithm", "solveMASP", "failed to solve master problem", 0);
			return 1;
		}
	}//END while loop

	writeStat(prob[0], cell, 0.0, probName);

	/*evaluating the optimal solution*/
	if (config.EVAL_FLAG == 1) {
		evaluate(stoc, prob, cell);
	}

	mem_free(observ);

	return 0;
}//END solveCell()

void writeStat(probType *prob, cellType *cell, double totRunTime, string probName) {
	FILE    *Stat;
	FILE    *StatTime;

	StatTime = openFile(outputDir, "time.dat", "a");
	fprintf(StatTime, "\nTotal (%d):\t", cell->k-1);

	fprintf(StatTime, "\ntotRunTime = %lf", totRunTime);

	fclose(StatTime);

	Stat = openFile(outputDir, "statistic.dat", "w");
	fprintf(Stat, "-------- Result summary for %s --------\n", probName);

	fprintf(Stat, "Samples simulated using SMPS stoch file.\n");
	fprintf(Stat, "A single aggregated cut added to master problem.\n\n");

	fprintf(Stat, "\nTotal number of iteration: %d\n", cell->k-1);
	fprintf(Stat, "Aggregated estimate at the incumbent solution = %lf\n", cell->incumbEst);
	fprintf(Stat, "\nIncumbent solution = ");
	printVector(cell->incumbX, prob->num->cols, Stat);

	fclose(Stat);
}//END WriteStat
