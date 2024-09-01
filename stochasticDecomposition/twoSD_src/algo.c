#include "twoSD.h"
#include <utils.h>

extern cString outputDir;
extern configType config;

int algo(oneProblem *orig, timeType *tim, stocType *stoc, cString inputDir, cString probName)
{
	dVector meanSol = NULL;
	probType **prob = NULL;
	cellType *cell1 = NULL;
	cellType *cell2 = NULL;
	batchSummary *batch = NULL;
	FILE *sFile1 = NULL, *sFile2 = NULL, *iFile = NULL, *bFile = NULL;
	dVector observ1 = NULL, observ2 = NULL;

	// Complete necessary initialization for the algorithm
	if (setupAlgo(orig, stoc, tim, &prob, &cell1, &batch, &meanSol))
	{
		goto TERMINATE;
	}

	// Duplicate cell1 to cell2
	cell2 = newCell(stoc, prob, meanSol);
	if (cell2 == NULL)
	{
		errMsg("allocation", "algo", "cell2", 0);
		goto TERMINATE;
	}

	printf("Starting two-stage stochastic decomposition.\n");

	char outputFileName1[100], outputFileName2[100];
	sprintf(outputFileName1, "%s_cell1.csv", probName);
	sprintf(outputFileName2, "%s_cell2.csv", probName);
	sFile1 = openFile(outputDir, outputFileName1, "w");
	sFile2 = openFile(outputDir, outputFileName2, "w");

	iFile = openFile(outputDir, "incumb.dat", "w");
	bFile = openFile(outputDir, "summary.dat", "w");

	printDecomposeSummary(bFile, probName, tim, prob);
	printDecomposeSummary(stdout, probName, tim, prob);

	// Allocate memory for observations
	if (!(observ1 = (dVector)arr_alloc(stoc->numOmega + 1, double)) ||
		!(observ2 = (dVector)arr_alloc(stoc->numOmega + 1, double)))
	{
		errMsg("allocation", "algo", "observ1 or observ2", 0);
		goto TERMINATE;
	}

	// Single while loop to process both cells
	while ((cell1->optFlag == false && cell1->k < config.MAX_ITER) ||
		   (cell2->optFlag == false && cell2->k < config.MAX_ITER))
	{
		// Generate antithetic samples for both cells using the same seed
		generateAntitheticOmegas(stoc, observ1, observ2, config.RUN_SEED);

		// Solve Cell 1 only if it hasn't reached optimality
		if (cell1->optFlag == false && solveCell(stoc, prob, cell1, observ1))
		{
			errMsg("algorithm", "algo", "failed to solve cell1 using 2-SD algorithm", 0);
			goto TERMINATE;
		}

		// Solve Cell 2 only if it hasn't reached optimality
		if (cell2->optFlag == false && solveCell(stoc, prob, cell2, observ2))
		{
			errMsg("algorithm", "algo", "failed to solve cell2 using 2-SD algorithm", 0);
			goto TERMINATE;
		}
	}

	// After both cells have been processed, print results

	// Print and write the results for Cell 1
	printf("\n\n\n\n Cell 1 Results...\n");
	printOptimizationSummary(cell1);
	writeOptimizationStatistics(sFile1, iFile, prob, cell1, 0);

	if (config.EVAL_FLAG == 1)
	{
		evaluate(sFile1, stoc, prob, cell1->subprob, cell1->incumbX);
	}
	else
	{
		fprintf(sFile1, "\n");
	}

	if (config.COMPROMISE_PROB)
	{
		buildCompromise(prob[0], cell1, batch);
	}

	// Print and write the results for Cell 2
	printf("\n\n\n\n Cell 2 Results...\n");
	printOptimizationSummary(cell2);
	writeOptimizationStatistics(sFile2, iFile, prob, cell2, 0);

	if (config.EVAL_FLAG == 1)
	{
		evaluate(sFile2, stoc, prob, cell2->subprob, cell2->incumbX);
	}
	else
	{
		fprintf(sFile2, "\n");
	}

	if (config.COMPROMISE_PROB)
	{
		buildCompromise(prob[0], cell2, batch);
	}

	if (config.COMPROMISE_PROB)
	{
		if (solveCompromise(prob[0], batch))
		{
			errMsg("algorithm", "algo", "failed to solve the compromise problem", 0);
			goto TERMINATE;
		}

		fprintf(bFile, "\n====================================================================================================================================\n");
		fprintf(bFile, "\n----------------------------------------- Compromise solution --------------------------------------\n\n");
		fprintf(stdout, "\n====================================================================================================================================\n");
		fprintf(stdout, "\n----------------------------------------- Compromise solution --------------------------------------\n\n");
		evaluate(sFile1, stoc, prob, cell1->subprob, batch->compromiseX);
		evaluate(sFile2, stoc, prob, cell2->subprob, batch->compromiseX);

		fprintf(bFile, "\n------------------------------------------- Average solution ---------------------------------------\n\n");
		fprintf(stdout, "\n------------------------------------------- Average solution ---------------------------------------\n\n");
		evaluate(sFile1, stoc, prob, cell1->subprob, batch->avgX);
		evaluate(sFile2, stoc, prob, cell2->subprob, batch->avgX);
	}

	fclose(sFile1);
	fclose(sFile2);
	fclose(iFile);
	fclose(bFile);
	printf("\nSuccessfully completed two-stage stochastic decomposition algorithm.\n");

	if (meanSol)
	{
		mem_free(meanSol);
	}
	freeBatchType(batch);
	freeCellType(cell1);
	freeCellType(cell2);
	freeProbType(prob, 2);
	mem_free(observ1);
	mem_free(observ2);
	return 0;

TERMINATE:
	if (meanSol)
	{
		mem_free(meanSol);
	}
	freeBatchType(batch);
	freeCellType(cell1);
	freeCellType(cell2);
	freeProbType(prob, 2);
	mem_free(observ1);
	mem_free(observ2);
	return 1;
}

int solveCell(stocType *stoc, probType **prob, cellType *cell, dVector observ)
{
	int m, omegaIdx, candidCut;
	bool newOmegaFlag;
	clock_t tic;

	// Main Algorithm
	tic = clock();

	// Step 1: Optimality tests
	if (optimal(prob, cell))
		return 0;

	cell->k++;

	// Step 2.c: Update omegaType with the latest observation
	omegaIdx = calcOmega(observ, 0, prob[1]->num->numRV, cell->omega, &newOmegaFlag, config.TOLERANCE);

	// Step 3: Solve the subproblem with candidate solution, form and update the candidate cut
	if ((candidCut = formSDCut(prob, cell, cell->candidX, omegaIdx, &newOmegaFlag, prob[0]->lb, CANDIDATE)) < 0)
	{
		errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
		return 1;
	}

	// Step 4: Solve subproblem with incumbent solution, and form an incumbent cut
	if (((cell->k - cell->iCutUpdt) % config.TAU == 0))
	{
		if ((formSDCut(prob, cell, cell->incumbX, omegaIdx, &newOmegaFlag, prob[0]->lb, INCUMBENT)) < 0)
		{
			errMsg("algorithm", "solveCell", "failed to create the incumbent cut", 0);
			return 1;
		}
	}

	// Step 5: Check improvement in predicted values at candidate solution
	if (!(cell->incumbChg) && cell->k > 1)
	{
		checkImprovement(prob[0], cell, candidCut);
	}

	// Step 6: Solve the master problem to obtain the new candidate solution
	if (solveQPMaster(prob[0]->num, prob[0]->dBar, cell, prob[0]->lb))
	{
		errMsg("algorithm", "solveCell", "failed to solve master problem", 0);
		return 1;
	}

	// Accumulate time
	cell->time.masterAccumTime += cell->time.masterIter;
	cell->time.subprobAccumTime += cell->time.subprobIter;
	cell->time.argmaxAccumTime += cell->time.argmaxIter;
	cell->time.optTestAccumTime += cell->time.optTestIter;
	cell->time.masterIter = cell->time.subprobIter = cell->time.optTestIter = cell->time.argmaxIter = 0.0;
	cell->time.iterTime = ((double)clock() - tic) / CLOCKS_PER_SEC;
	cell->time.iterAccumTime += cell->time.iterTime;

	return 0;
}
