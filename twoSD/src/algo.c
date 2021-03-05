/*
 * algo.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *
 *
 *   Edited on: Dec 14, 2020 as part of the SD-integer project
 *      Author: Siavash Tabrizian
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to stabrizian (at) smu (dot) edu
 *
 */

#include "bnc.h" 

extern cString outputDir;
extern configType config;


int algo(oneProblem *orig, timeType *tim, stocType *stoc, cString inputDir, cString probName) {
	dVector	 meanSol = NULL;
	probType **prob = NULL;
	cellType *cell = NULL;
	batchSummary *batch = NULL;
	FILE 	*sFile = NULL, *iFile = NULL, *dFile = NULL, *oFile = NULL;

	/* open solver environment */
	openSolver();

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell, &batch, &meanSol, PROB_QP) )
		goto TERMINATE;

#if defined(LPMIP_PRINT)
	writeProblem(orig->lp, "meanvalueprob.lp");
#endif


	printf("Starting two-stage stochastic decomposition.\n");
	dFile = openFile(outputDir, "detailedResults.csv", "w");
	sFile = openFile(outputDir, "results.dat", "w");
	iFile = openFile(outputDir, "incumb.dat", "w");
	oFile = openFile(outputDir, "omega.dat", "w");
	printDecomposeSummary(sFile, probName, tim, prob);
	printDecomposeSummary(stdout, probName, tim, prob);

	for ( int rep = 0; rep < config.MULTIPLE_REP; rep++ ) {
		fprintf(sFile, "\n====================================================================================================================================\n");
		fprintf(sFile, "Replication-%d\n", rep+1);
		fprintf(stdout, "\n====================================================================================================================================\n");
		fprintf(stdout, "Replication-%d\n", rep+1);

		/* setup the seed to be used in the current iteration */
		config.RUN_SEED[0] = config.RUN_SEED[rep+1];
		config.EVAL_SEED[0] = config.EVAL_SEED[rep+1];

		if ( rep != 0 ) {
			/* clean up the cell for the next replication */
			if ( cleanCellType(cell, prob[0], meanSol) ) {
				errMsg("algorithm", "algo", "failed clean the problem cell", 0);
				goto TERMINATE;
			}
		}

		clock_t tic = clock();

#if  defined(experiment)
		/* evaluating a test solution*/
		double eval2;
		double candidX[11] = { 30.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0 };
		evaluate(sFile, stoc, prob, cell->subprob, candidX);
		goto TERMINATE;
#endif //  defined(experiment)

		/* Use two-stage stochastic decomposition algorithm to solve the problem */
		if ( config.MASTER_TYPE == PROB_MILP ) {
			if (branchbound(stoc, prob, cell, -INFINITY, INFINITY)) {
				errMsg("algorithm", "algo", "failed to solve the problem using integer B&C 2-SD algorithm", 0);
				goto TERMINATE;
			}
		}
		else {
			if ( solveCell(stoc, prob, cell) ) {
				errMsg("algorithm", "algo", "failed to solve the cell using 2-SD algorithm", 0);
				goto TERMINATE;
			}
		}

		cell->time.repTime = ((double) clock() - tic)/CLOCKS_PER_SEC;

		/* Write solution statistics for optimization process */
		if (rep == 0 ) {
			writeOptimizationSummary(sFile, iFile, prob, cell, true);
			writeOptimizationSummary(stdout, NULL, prob, cell, true);
		}
		else {
			writeOptimizationSummary(sFile, iFile, prob, cell, false);
			writeOptimizationSummary(stdout, NULL, prob, cell, false);
		}
		writeOptimizationStatistics(dFile, iFile, prob, cell, rep);

		/* evaluate the optimal solution*/
		if (config.EVAL_FLAG == 1)
		{
			//LB estimate: CTx (dot product) + height of the cut(subroutine: )
			evaluate(dFile, sFile, stoc, prob, cell->subprob, cell->incumbX);
		}
		else
			fprintf(dFile, "\n");
			
		/* print the omega structure */
		writeOmegaHist(oFile, cell, rep);

		/* Save the batch details and build the compromise problem. */
		if ( config.COMPROMISE_PROB) {
			buildCompromise(prob[0], cell, batch);
		}
	}

	if ( config.COMPROMISE_PROB ) {
		/* Solve the compromise problem. */
		if ( solveCompromise(prob[0], batch)) {
			errMsg("algorithm", "algo", "failed to solve the compromise problem", 0);
			goto TERMINATE;
		}

		fprintf(sFile, "\n====================================================================================================================================\n");
		fprintf(sFile, "\n----------------------------------------- Compromise solution --------------------------------------\n\n");
		fprintf(sFile, "\n====================================================================================================================================\n");
		fprintf(sFile, "\n----------------------------------------- Compromise solution --------------------------------------\n\n");
		/* Evaluate the compromise solution */
		evaluate(NULL, sFile, stoc, prob, cell->subprob, batch->compromiseX);

		fprintf(sFile, "\n------------------------------------------- Average solution ---------------------------------------\n\n");
		fprintf(stdout, "\n------------------------------------------- Average solution ---------------------------------------\n\n");
		/* Evaluate the average solution */
		evaluate(NULL,NULL, stoc, prob, cell->subprob, batch->avgX);
	}

	fclose(sFile); fclose(iFile); fclose(dFile);
	printf("\nSuccessfully completed two-stage stochastic decomposition algorithm.\n");

	/* free up memory before leaving */
	if (meanSol) mem_free(meanSol);
	freeBatchType(batch);
	freeCellType(cell);
	freeProbType(prob, 2);
	return 0;

	TERMINATE:
	if(meanSol) mem_free(meanSol);
	freeBatchType(batch);
	freeCellType(cell);
	freeProbType(prob, 2);
	return 1;
}//END algo()

void writeOptimizationSummary(FILE *soln, FILE *incumb, probType **prob, cellType *cell, bool header) {

	if ( header ) {
		fprintf(soln, "\n--------------------------------------- Problem Information ----------------------------------------\n\n");
		fprintf(soln, "Problem                                : %s\n", prob[0]->name);
		fprintf(soln, "First Stage Rows                       : %d\n", prob[0]->num->rows);
		fprintf(soln, "First Stage Columns                    : %d\n", prob[0]->num->cols);
	}

	fprintf(soln, "\n------------------------------------------- Optimization -------------------------------------------\n\n");

	fprintf(soln, "Algorithm                              : Two-stage Stochastic Decomposition\n");
	fprintf(soln, "Number of iterations                   : %d\n", cell->k);
	fprintf(soln, "Lower bound estimate                   : %f\n", cell->incumbEst);
	fprintf(soln, "Total time                             : %f\n", cell->time.repTime);
	fprintf(soln, "Total time to solve master             : %f\n", cell->time.masterAccumTime);
	fprintf(soln, "Total time to solve subproblems        : %f\n", cell->time.subprobAccumTime);
	fprintf(soln, "Total time in argmax procedure         : %f\n", cell->time.argmaxAccumTime);
	fprintf(soln, "Total time in verifying optimality     : %f\n", cell->time.optTestAccumTime);

	if ( incumb != NULL ) {
		printVector(cell->incumbX, prob[0]->num->cols, incumb);
	}

}//END WriteStat
