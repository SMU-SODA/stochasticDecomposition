/*
 * evaluate.c
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

/* Haven't intergrated the multi-cut and aggregated-cut versions */
void evaluate(stocType *stoc, probType **prob, cellType *cell) {
//	vector 	observ;
//	double 	obj, mean, variance, stdev, temp, CI[2], totalTime;
//	int		status, cnt, stat2, i, m, sOffset = 0;
//	//    clock_t	tic, toc;
//	FILE	*ePtr;
//
//	if ( !(observ = (vector) arr_alloc(stoc->numOmega + 1, double)) )
//		errMsg("allocation", "evaluateOpt", "observ", 0);
//
//	printf("\n\nEvaluating optimal solution");
//	//tic = clock();
//
//
//	/* open a file to record evaluation results */
//	ePtr = openFile(outputDir, "eval.dat", "w");
//
//	for (i = 1; i <numAgents; i++) {
//		cnt = 0.0; mean = 0.0; variance = 0.0; stdev = INFBOUND; cnt = 0;
//		chgRHSwMean(prob[i]->bBar, prob[i]->Cbar, cell->spRHS, cell->incumbX);
//
//		while (cnt < config.EVAL_MIN_ITER ) {
//			/* use the stoc file to generate observations */
//			generateOmega(stoc, observ, &config.EVAL_SEED);
//
//			for ( m = 0; m < stoc->numOmega; m++ )
//				observ[m] -= stoc->mean[m];          /* store the mean rv in observ */
//
//			/* setup and solve subproblem */
//			status = chgRHSwRand(cell->subprob->lp, prob[i]->num, prob[i]->coord, observ + sOffset -1, cell->spRHS, cell->incumbX);
//			if ( status ) {
//				errMsg("algorithm", "evaluateOpt", "failed to setup the subproblem",0);
//				exit(1);
//			}
//
//#if 0
//			int     status;
//			char probName[NAMESIZE];
//			sprintf(probName,"evaluation_%d.lp", cell[i]->ID);
//			status = writeProblem(cell[i]->sp->lp, probName);
//			if ( status ) {
//				errMsg("write problem", "new_subprob", "failed to write subproblems problem to file",0);
//				exit(1);
//			}
//#endif
//			status = solveProblem(cell->subprob->lp, cell->subprob->name, cell->subprob->type, &stat2);
//			if ( status ) {
//				if ( stat2 == STAT_INFEASIBLE ) {
//					/* subproblem is infeasible */
//					printf("Warning:: Subproblem is infeasible: need to create feasibility cut.\n");
//					exit(1);
//				}
//				else {
//					errMsg("algorithm", "evaluateOpt", "failed to solve subproblem in solver", 0);
//					exit(1);
//				}
//			}
//
//			/* use subproblem objective and compute evaluation statistics */
//			obj = getObjective(cell->subprob->lp, PROB_LP);
//
//			if ( cnt == 0 )
//				mean = obj;
//			else {
//				temp = mean;
//				mean = mean + (obj - mean) / (double) (cnt + 1);
//				variance  = (1 - 1 / (double) cnt) * variance
//						+ (cnt + 1) * (mean - temp) * (mean - temp);
//				stdev = sqrt(variance/ (double) cnt);
//			}
//			cnt++;
//			/* Print the results every once in a while for long runs */
//			if (!(cnt % 100)) {
//				printf(".");
//				fflush(stdout);
//			}
//			if (!(cnt % 10000))
//				printf("\n\nObs:%d mean:%lf   error: %lf \n 0.90 CI: [%lf , %lf]\n", cnt, mean, 3.29 * stdev / mean,
//						mean - 1.645 * stdev, mean + 1.645 * stdev);
//		}//END while loop
//		CI[0] = mean - 1.645 * stdev;
//		CI[1] = mean + 1.645 * stdev;
//
//		/* Print the value of the solution to a file and the screen */
//		printf("\nFinal Estimate :: obs:%d, mean:%lf, 0.90 C.I.: [%lf , %lf] \n\n",cnt, mean, CI[0], CI[1]);
//
//		fprintf(ePtr, "\t\t%lf\t%lf\t[%lf, %lf]\ttime:%lf\t%d\n",mean, stdev, CI[0], CI[1], totalTime, cnt);
//
//		if (!(prob[i]->omegas == NULL))
//			sOffset += prob[i]->omegas->numRV;
//
//	}//END agents for loop
//
//	fclose(ePtr);
//	mem_free(observ);
}
