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
		fprintf(soln, "Replication,Iterations,LB estimate,Total time,Master time,Subproblem time,Optimality time,Argmax time,Reduce time,"
			"depth,tot nodes,d nodes, i nodes, maxiter,UB Estimate,Error,CI-L,CI-U,outcome\n");

	fprintf(soln, "%d,%d,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%d,%d,%d,%d,%d", rep+1, cell->k, cell->incumbEst,cell->time.repTime, cell->time.masterAccumTime,
			cell->time.subprobAccumTime, cell->time.optTestAccumTime, cell->time.argmaxAccumTime, 0.0, cell->depth,cell->tot_nodes,cell->d_nodes,cell->int_nodes,cell->maxiter_nodes);

	printVector(cell->incumbX, prob[0]->num->cols, incumb);

}//END writeOptimizationStatistics()

void writeEvaluationStatistics(FILE *soln, double mean, double stdev, int cnt) {

	fprintf(soln, ",%.4lf,%.4lf,%.4lf,%.4lf,%d\n", mean, 3.29 * stdev / mean, mean - 1.645 * stdev, mean + 1.645 * stdev, cnt);

}//END writeEvaluationSummary()


void writeOmegaHist(FILE *soln, cellType *cell, int rep) {

	/* Print header for the first replication*/
	if (rep == 0)
	{
		fprintf(soln, "Frequency table:\n");
		fprintf(soln, "Sen: %-5s"," ");
		for (int o = 0; o < cell->omega->cnt; o++)
		{
			fprintf(soln, "\t%-5d", o+1);
		}
		fprintf(soln, "\n     %-5s", "-");
		for (int o = 0; o < cell->omega->cnt; o++)
		{
			fprintf(soln, "\t%-5s", "---");
		}
	}
	fprintf(soln, "\nrep: %-5d|",rep+1);
	for (int o = 0; o < cell->omega->cnt; o++)
	{
		fprintf(soln, "\t%-5d", cell->omega->weights[o]);
	}

}//END writeOmegaHist()

/* Prints summary statistics to the console output. */
//void printOptimizationSummary(cellType *cell) {
//
//	fprintf(stdout, "\n------------------------------------------------------------ Optimization ---------------------------------------------------------\n");
//	if ( config.MASTER_TYPE == PROB_QP )
//		fprintf(stdout, "Algorithm                          : Regularized Benders Decomposition\n");
//	else
//		fprintf(stdout, "Algorithm                          : Benders Decomposition\n");
//	fprintf(stdout, "Number of iterations               : %d", cell->k);
//	if ( cell->k == config.MAX_ITER)
//		fprintf(stdout, "*\n");
//	else
//		fprintf(stdout, "\n");
//	fprintf(stdout, "Lower bound estimate               : %f\n", cell->incumbEst);
//	fprintf(stdout, "Total time                         : %f\n", cell->time->repTime);
//	fprintf(stdout, "Total time to solve master         : %f\n", cell->time->masterAccumTime);
//	fprintf(stdout, "Total time to solve subproblems    : %f\n", cell->time->subprobAccumTime);
//	fprintf(stdout, "Total time to verify optimality    : %f\n", cell->time->optTestAccumTime);
//	if ( config.SAMPLED_SP ) {
//		fprintf(stdout, "Total time for argmax operation    : %f\n", cell->time->argmaxAccumTime);
//	}
//	if ( config.SCENRED ) {
//		fprintf(stdout, "Total time for scenario reduction  : %f\n", cell->time->reduceTime);
//	}
//
//}//END writeOptimizationSummary()
//
//void printEvaluationSummary(FILE *soln, double mean, double stdev, int cnt) {
//
//	/* Write the evaluation results to the summary file */
//	fprintf(soln, "\n-------------------------------------------- Evaluation --------------------------------------------\n\n");
//	fprintf(soln, "Upper bound estimate               : %lf\n", mean);
//	fprintf(soln, "Error in estimation                : %lf\n", 3.29 * stdev / mean);
//	fprintf(soln, "Confidence interval at 95%%         : [%lf, %lf]\n", mean - 1.645 * stdev, mean + 1.645 * stdev);
//	fprintf(soln, "Number of observations             : %d\n", cnt);
//
//
//}//END writeEvaluationSummary()

