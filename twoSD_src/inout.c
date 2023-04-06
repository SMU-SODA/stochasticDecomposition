/*
 * inout.c
 *
 *  Created on: Feb 13, 2019
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"

extern configType config;

void writeOptimizationStatistics(FILE *soln, FILE *incumb, probType **prob, cellType *cell, int rep) {

	/* Print header for the first replication*/
	if ( rep == 0)
		fprintf(soln, "Replication\tIterations\tLB estimate\tTotal time\tMaster time\t Subproblem time\t Optimality time\tArgmax time\t"
				"UB Estimate\tError\tCI-L\tCI-U\tOutcomes\n");

	fprintf(soln, "%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf", rep+1, cell->k, cell->incumbEst, cell->time.repTime,
			cell->time.masterAccumTime, cell->time.subprobAccumTime, cell->time.optTestAccumTime, cell->time.argmaxAccumTime);

	if ( incumb != NULL ) {
		if ( config.MASTER_TYPE == PROB_QP )
			printVector(cell->incumbX, prob[0]->num->cols, incumb);
		else
			printVector(cell->candidX, prob[0]->num->cols, incumb);
	}

}//END writeOptimizationStatistics()

void writeEvaluationStatistics(FILE *soln, double mean, double stdev, int cnt) {

	fprintf(soln, "\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%d\n", mean, 3.29 * stdev / mean, mean - 1.645 * stdev, mean + 1.645 * stdev, cnt);

}//END writeEvaluationSummary()

/* Prints summary statistics to the console output. */
void printOptimizationSummary(cellType *cell) {

	fprintf(stdout, "\n------------------------------------------------------------ Optimization ---------------------------------------------------------\n");
	fprintf(stdout, "Algorithm                          : Two-stage Stochastic Decomposition\n");
	fprintf(stdout, "Number of iterations               : %d", cell->k);
	if ( cell->k == config.MAX_ITER)
		fprintf(stdout, "*\n");
	else
		fprintf(stdout, "\n");
	fprintf(stdout, "Number of unique observations      : %d\n", cell->omega->cnt);
	fprintf(stdout, "Lower bound estimate               : %f\n", cell->incumbEst);
	fprintf(stdout, "Total time                         : %f\n", cell->time.repTime);
	fprintf(stdout, "Total time to solve master         : %f\n", cell->time.masterAccumTime);
	fprintf(stdout, "Total time to solve subproblems    : %f\n", cell->time.subprobAccumTime);
	fprintf(stdout, "Total time to verify optimality    : %f\n", cell->time.optTestAccumTime);
	fprintf(stdout, "Total time for argmax operation    : %f\n", cell->time.argmaxAccumTime);

}//END writeOptimizationSummary()

void printEvaluationSummary(FILE *soln, double mean, double stdev, int cnt) {

	/* Write the evaluation results to the summary file */
	fprintf(soln, "\n-------------------------------------------- Evaluation --------------------------------------------\n\n");
	fprintf(soln, "Upper bound estimate               : %lf\n", mean);
	fprintf(soln, "Error in estimation                : %lf\n", 3.29 * stdev / mean);
	fprintf(soln, "Confidence interval at 95%%         : [%lf, %lf]\n", mean - 1.645 * stdev, mean + 1.645 * stdev);
	fprintf(soln, "Number of observations             : %d\n", cnt);


}//END writeEvaluationSummary()

