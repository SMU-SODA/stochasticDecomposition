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

extern cString outputDir;
extern configType config;

dVector roundX(dVector a, int len);
int isCandidInt(dVector candidate, int size);

int algo(oneProblem *orig, timeType *tim, stocType *stoc, cString inputDir, cString probName) {
	dVector	 meanSol = NULL;
	probType **prob = NULL;
	cellType *cell = NULL;
	probType **clone_prob = NULL;
	cellType *clone_cell = NULL;
	dVector	lb = NULL;
	batchSummary *batch = NULL;
	FILE 	*sFile = NULL, *iFile = NULL;

	/* open solver environment */
	openSolver();

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell, &batch, &meanSol, &lb) )
		goto TERMINATE;

#if defined(LPMIP_PRINT)
	writeProblem(orig->lp, "meanvalueprob.lp");
#endif

	/* create clone problems */
	if (setupClone(orig, stoc, tim, &clone_prob, &clone_cell, &meanSol, &lb))
		goto TERMINATE;

	printf("Starting two-stage stochastic decomposition.\n");
	sFile = openFile(outputDir, "results.dat", "w");
	iFile = openFile(outputDir, "incumb.dat", "w");
	printDecomposeSummary(sFile, probName, tim, prob);
	printDecomposeSummary(stdout, probName, tim, prob);

	for ( int rep = 0; rep < config.NUM_REPS; rep++ ) {
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
		/* Use two-stage stochastic decomposition algorithm to solve the problem */
		if ( solveCell(stoc, prob, cell,clone_prob,clone_cell) ) {
			errMsg("algorithm", "algo", "failed to solve the cell using 2-SD algorithm", 0);
			goto TERMINATE;
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

		/* evaluate the optimal solution*/
		if (config.EVAL_FLAG == 1)
		{
			//LB estimate: CTx (dot product) + height of the cut(subroutine: )
			//cell->incumbEst = vXvSparse(cell->incumbX, prob[0]->dBar) + maxCutHeight(cell->cuts, cell->sampleSize, cell->incumbX, prob[0]->num->cols, cell->lb);
			evaluate(sFile, stoc, prob, cell->subprob, cell->incumbX);
		}
			

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
		evaluate(sFile, stoc, prob, cell->subprob, batch->compromiseX);

		fprintf(sFile, "\n------------------------------------------- Average solution ---------------------------------------\n\n");
		fprintf(stdout, "\n------------------------------------------- Average solution ---------------------------------------\n\n");
		/* Evaluate the average solution */
		evaluate(sFile, stoc, prob, cell->subprob, batch->avgX);
	}

	fclose(sFile); fclose(iFile);
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
	freeCellType(clone_cell);
	freeCellType(cell);
	freeProbType(clone_prob, 2);
	freeProbType(prob, 2);
	mem_free(lb);
	return 1;
}//END algo()

int intalgo(oneProblem *orig, timeType *tim, stocType *stoc, cString inputDir, cString probName) {
	dVector	 meanSol = NULL;
	probType **prob = NULL;
	cellType *cell = NULL;
	probType **clone_prob = NULL;
	cellType *clone_cell = NULL;
	dVector	lb = NULL;
	batchSummary *batch = NULL;
	FILE 	*sFile = NULL, *iFile = NULL;

	/* open solver environment */
	openSolver();

	/* complete necessary initialization for the algorithm */
	if (setupAlgo(orig, stoc, tim, &prob, &cell, &batch, &meanSol, lb))
		goto TERMINATE;

	/* create clone problems */
	if (setupClone(orig, stoc, tim, &clone_prob, &clone_cell, &batch, &meanSol,&lb))
		goto TERMINATE;

	printf("Starting two-stage stochastic decomposition.\n");
	sFile = openFile(outputDir, "results.dat", "w");
	iFile = openFile(outputDir, "incumb.dat", "w");
	printDecomposeSummary(sFile, probName, tim, prob);
	printDecomposeSummary(stdout, probName, tim, prob);

	for (int rep = 0; rep < config.NUM_REPS; rep++) {
		fprintf(sFile, "\n====================================================================================================================================\n");
		fprintf(sFile, "Replication-%d\n", rep + 1);
		fprintf(stdout, "\n====================================================================================================================================\n");
		fprintf(stdout, "Replication-%d\n", rep + 1);

		/* setup the seed to be used in the current iteration */
		config.RUN_SEED[0] = config.RUN_SEED[rep + 1];
		config.EVAL_SEED[0] = config.EVAL_SEED[rep + 1];

		if (rep != 0) {
			/* clean up the cell for the next replication */
			if (cleanCellType(cell, prob[0], meanSol)) {
				errMsg("algorithm", "algo", "failed clean the problem cell", 0);
				goto TERMINATE;
			}
		}

		if (changeQPSolverType(1))
		{
			errMsg("algorithm", "algo", "failed to set primal algorithm for QP", 0);
			goto TERMINATE;
		}

		clock_t tic = clock();
		/* Use two-stage stochastic decomposition algorithm to solve the problem */
		if (solveCell(stoc, prob, cell, clone_prob, clone_cell)) {
			errMsg("algorithm", "algo", "failed to solve the cell using 2-SD algorithm", 0);
			goto TERMINATE;
		}
		printf("\nProblem type before: %i \n", getProbType(cell->master->lp));
		printf("objective before: %4.6lf \n", getObjective(cell->master->lp, 5));
		//set sigma = 0 which is in master.c
		//changerihgthand side for changing the alpha incumbents master.c (changeQPRHS)

		printf("\n");
		printVector(cell->candidX, prob[0]->num->cols, stdout);

		fprintf(sFile, "\nAdding GMI and MIR cuts\n\n");
		fprintf(stdout, "\nAdding GMI and MIR cuts \n\n");

		/* Use GMI and MIR cutting planes to solve the SD-optimized problem */
		if (solveIntCell(stoc, prob, cell)) {
			errMsg("algorithm", "algo", "failed to solve the cell using GMI and MIR algorithm", 0);
			goto TERMINATE;
		}
		fprintf(sFile, "\n %i GMI and %i MIR cuts are added \n\n", cell->GMIcuts->cnt, cell->MIRcuts->cnt);
		fprintf(stdout, "\n %i GMI and %i MIR are added \n\n",cell->GMIcuts->cnt, cell->MIRcuts->cnt);
		printf("\n");
		printVector(cell->candidX, prob[0]->num->cols, stdout);

		cell->time.repTime = ((double)clock() - tic) / CLOCKS_PER_SEC;

		/* Write solution statistics for optimization process */
		if (rep == 0) {
			writeOptimizationSummary(sFile, iFile, prob, cell, true);
			writeOptimizationSummary(stdout, NULL, prob, cell, true);
		}
		else {
			writeOptimizationSummary(sFile, iFile, prob, cell, false);
			writeOptimizationSummary(stdout, NULL, prob, cell, false);
		}

		/* evaluate the optimal solution*/
		if (config.EVAL_FLAG == 1)
			evaluate(sFile, stoc, prob, cell->subprob, cell->incumbX);

		/* Save the batch details and build the compromise problem. */
		if (config.COMPROMISE_PROB) {
			buildCompromise(prob[0], cell, batch);
		}
	}

	if (config.COMPROMISE_PROB) {
		/* Solve the compromise problem. */
		if (solveCompromise(prob[0], batch)) {
			errMsg("algorithm", "algo", "failed to solve the compromise problem", 0);
			goto TERMINATE;
		}

		fprintf(sFile, "\n====================================================================================================================================\n");
		fprintf(sFile, "\n----------------------------------------- Compromise solution --------------------------------------\n\n");
		fprintf(sFile, "\n====================================================================================================================================\n");
		fprintf(sFile, "\n----------------------------------------- Compromise solution --------------------------------------\n\n");
		/* Evaluate the compromise solution */
		evaluate(sFile, stoc, prob, cell->subprob, batch->compromiseX);

		fprintf(sFile, "\n------------------------------------------- Average solution ---------------------------------------\n\n");
		fprintf(stdout, "\n------------------------------------------- Average solution ---------------------------------------\n\n");
		/* Evaluate the average solution */
		evaluate(sFile, stoc, prob, cell->subprob, batch->avgX);
	}

	fclose(sFile); fclose(iFile);
	printf("\nSuccessfully completed two-stage stochastic decomposition algorithm.\n");

	/* free up memory before leaving */
	if (meanSol) mem_free(meanSol);
	freeBatchType(batch);
	freeCellType(cell);
	freeProbType(prob, 2);
	return 0;

TERMINATE:
	if (meanSol) mem_free(meanSol);
	freeBatchType(batch);
	freeCellType(cell);
	freeProbType(prob, 2);
	mem_free(lb);
	return 1;
}//END algo()

int solveCell(stocType *stoc, probType **prob, cellType *cell, probType **clone_prob, cellType *clone_cell) {
	dVector 	observ;
	clock_t		tic;
	bool breakLoop = false;



	/* -+-+-+-+-+-+-+-+-+-+-+-+-+-+- Main Algorithm -+-+-+-+-+-+-+-+-+-+-+-+-+-+- */
	observ = (dVector) arr_alloc(stoc->numOmega + 1, double);

	/******* 0. Initialization: The algorithm begins by solving the master problem as a QP *******/
	while (cell->optFlag == false && cell->k < config.MAX_ITER) {
		
		tic = clock();

		if (mainloopSDCell(stoc, prob, cell, &breakLoop, observ)) {
			errMsg("Callback", "usersolve", "failed to solve Benders cell for the node problem", 0);
			goto TERMINATE;
		}

		cell->time.masterAccumTime += cell->time.masterIter; cell->time.subprobAccumTime += cell->time.subprobIter;
		cell->time.argmaxAccumTime += cell->time.argmaxIter; cell->time.optTestAccumTime += cell->time.optTestIter;
		cell->time.masterIter = cell->time.subprobIter = cell->time.optTestIter = cell->time.argmaxIter = 0.0;
		cell->time.iterTime = ((double) clock() - tic)/CLOCKS_PER_SEC; cell->time.iterAccumTime += cell->time.iterTime;
		printf("*"); fflush(stdout);

		if (breakLoop)
			break;

	}//END while loop

#if defined(LPMIP_PRINT)
	 // Print before branch and bound -------------------------
	printf("\nQP relaxation solution:\n");
	printVector(cell->incumbX, cell->master->mac - 1, NULL);
	printf("\nQP relaxation Estimate %0.4f", cell->incumbEst);
	writeProblem(cell->master->lp, "finalMaster_QP.lp");
	//--------------------------------------------------------
#endif // defined(LPMIP_PRINT)


#if defined(PHASE1ANLYS)
	if (phase_one_analysis(stoc, prob, cell)) {
		errMsg("algorithm", "solveCell-phase-one-analysis", "failed to run analysis after the phase one", 0);
		goto TERMINATE;
	}
#endif // defined(PHASE1ANLYS)

	//2a - change master to MILP solve using callback 

	/* Clone the current cell to create an LP cell problem */

	/* Copy the cell to the cloned cell. */
	if (copyCell(cell,clone_cell)) {
		errMsg("algo", "copyCell", "failed to create a copy of the cell", 0);
		return 1;
	}
	clone_cell->incumbEst = prob[0]->lb;
	clone_cell->candidEst = prob[0]->lb;
	clone_cell->incumbX = roundX(clone_cell->candidX, prob[0]->num->cols);

	if (config.SMIP != MILP) {
		/* Turn the clone problem to LP */
		QPtoLP(stoc, prob, clone_cell, 1);

		/* Phase-1 has completed, we have an approximation obtained by solving the relaxed problem.
		* Phase-2 be used to impose integrality through costom procedures or the callback.  */
		clone_cell->optFlag = false;
		if (bendersCallback(stoc, prob, clone_cell)) {
			errMsg("algorithm", "solverCell", "failed to run the Benders-callback routine", 0);
			goto TERMINATE;
		}
		/* Update the incumbent estimate */
		clone_cell->incumbEst = vXvSparse(clone_cell->incumbX, prob[0]->dBar) + maxCutHeight(clone_cell->cuts, clone_cell->sampleSize, clone_cell->incumbX, prob[0]->num->cols, clone_cell->lb);

	}

#if defined(LPMIP_PRINT)
	printVector(clone_cell->incumbX, clone_cell->master->mac - 1, NULL);
	printf("\nEstimate %0.4f", clone_cell->incumbEst);
#endif // defined(LPMIP_PRINT)

	cell->incumbEst = clone_cell->incumbEst;
	*cell->incumbX = *clone_cell->incumbX;

	mem_free(observ);
	return 0;

	TERMINATE:
	mem_free(observ);
	return 1;
}//END solveCell()

int mainloopSDCell(stocType *stoc, probType **prob, cellType *cell, bool *breakLoop, dVector observ)
{
	int			m, candidCut, obs;
	
	cell->k++;
#if defined(STOCH_CHECK) || defined(ALGO_CHECK)
	printf("\nIteration-%d :: \n", cell->k);
#else
	if ((cell->k - 1) % 100 == 0) {
		printf("\nIteration-%4d: ", cell->k);
	}
#endif



	/******* 1. Optimality tests *******/
	if (optimal(prob, cell))
	{
		(*breakLoop) = true; return 0;
	}

	/******* 2. Generate new observations, and add it to the set of observations *******/
	cell->sampleSize += config.SAMPLE_INCREMENT;
	for (obs = 0; obs < config.SAMPLE_INCREMENT; obs++) {
		/* (a) Use the stoc file to generate observations */
		generateOmega(stoc, observ, config.TOLERANCE, &config.RUN_SEED[0], NULL);

		/* (b) Since the problem already has the mean values on the right-hand side, remove it from the original observation */
		for (m = 0; m < stoc->numOmega; m++)
			observ[m] -= stoc->mean[m];

		/* (d) update omegaType with the latest observation. If solving with incumbent then this update has already been processed. */
		cell->sample->omegaIdx[obs] = calcOmega(observ - 1, 0, prob[1]->num->numRV, cell->omega, &cell->sample->newOmegaFlag[obs], config.TOLERANCE);
	}

	/******* 3. Solve the subproblem with candidate solution, form and update the candidate cut *******/
	if ((candidCut = formSDCut(prob, cell, cell->candidX, prob[0]->lb, 0)) < 0) {
		errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
		return 1;
	}

	/******* 4. Solve subproblem with incumbent solution, and form an incumbent cut *******/
	if (((cell->k - cell->iCutUpdt) % config.TAU == 0)) {
		if ((cell->iCutIdx = formSDCut(prob, cell, cell->incumbX, prob[0]->lb, 0)) < 0) {
			errMsg("algorithm", "solveCell", "failed to create the incumbent cut", 0);
			return 1;
		}
		cell->iCutUpdt = cell->k;
	}

	/******* 5. Check improvement in predicted values at candidate solution *******/
	if (!(cell->incumbChg) && cell->k > 1)
		/* If the incumbent has not changed in the current iteration */
		checkImprovement(prob[0], cell, candidCut);

	/******* 6. Solve the master problem to obtain the new candidate solution */
	if (solveQPMaster(prob[0]->num, prob[0]->dBar, cell, prob[0]->lb)) {
		errMsg("algorithm", "solveCell", "failed to solve master problem", 0);
		return 1;
	}

	return 0;

}

int phase_one_analysis(stocType *stoc, probType **prob, cellType *cell)
{
	int status = 0;
	printf("\n--------------------------------Phase 1 summary------------------------------------------\n");
	//0a - (LB on xLP is cell->incumbEst)
	//0b - xLP is phase1 solution -> evaluate xLP (UB) 
	/* Evaluate the compromise solution */
	printf("LB estimate: %0.4f\n", cell->incumbEst);
	evaluate(NULL, stoc, prob, cell->subprob, cell->candidX);
	

	//1a - Change the master to MILP no callback solve to optimality - > xIP (LB)
	//1b - check the optimality of (xIP) using pretest 
	//1c - Evaluate the solution (xIP) using evaluate(sFile, stoc, prob, cell->subprob, cell->incumbX); -> (UB)
	/* Turn the clone problem to LP */
	QPtoLP(stoc, prob, cell, 0);
	/* Launch the solver to solve the MIP master problem */
	if (solveProblem(cell->master->lp, cell->master->name, cell->master->type, cell->master->mar, cell->master->mac,
		&status, config.SMIP_OPTGAP)) {
		errMsg("algorithm", "algo-after-phase1", "failed to solve the master problem", 0);
		return 1;
	}
	/* Get the most recent optimal solution to master program */
	if (getPrimal(cell->master->lp, cell->candidX, prob[0]->num->cols)) {
		errMsg("algorithm", "solveMaster", "failed to obtain the primal solution for master", 0);
		return 1;
	}

	/* Find the highest cut at the candidate solution. where cut_height = alpha - beta . x */ //
	cell->candidEst = vXvSparse(cell->candidX, prob[0]->dBar) + maxCutHeight(cell->cuts, cell->sampleSize, cell->candidX, prob[0]->num->cols, prob[0]->lb);
	printf("\n\nLP estimate: %0.4f\n", cell->candidEst);
	printVector(cell->candidX, prob[0]->num->cols, NULL);
	evaluate(NULL, stoc, prob, cell->subprob, cell->candidX);

	/* Turn the clone problem to MILP */
	LPtoMILP(stoc, prob, cell);
	/* Launch the solver to solve the MIP master problem */
	if (solveProblem(cell->master->lp, cell->master->name, cell->master->type, cell->master->mar, cell->master->mac,
		&status, config.SMIP_OPTGAP)) {
		errMsg("algorithm", "algo-after-phase1", "failed to solve the master problem", 0);
		return 1;
	}
	/* Get the most recent optimal solution to master program */
	if (getPrimal(cell->master->lp, cell->candidX, prob[0]->num->cols)) {
		errMsg("algorithm", "solveMaster", "failed to obtain the primal solution for master", 0);
		return 1;
	}

	/* Find the highest cut at the candidate solution. where cut_height = alpha - beta . x */ //
	cell->candidEst = vXvSparse(cell->candidX, prob[0]->dBar) + maxCutHeight(cell->cuts, cell->sampleSize, cell->candidX, prob[0]->num->cols, prob[0]->lb);
	printf("\n\nMILP estimate: %0.4f\n", cell->candidEst);
	printVector(cell->candidX, prob[0]->num->cols, NULL);
	evaluate(NULL, stoc, prob, cell->subprob, cell->candidX);

	//summary for the second one 

	printf("----------------------------------------------------------------------------------------------------\n");

	return 0;
}

int mainloopSDCell_callback(stocType *stoc, probType **prob, cellType *cell, bool *breakLoop, dVector observ)
{
	int			m, candidCut, obs;

	cell->k++;
#if defined(STOCH_CHECK) || defined(ALGO_CHECK)
	printf("\nIteration-%d :: \n", cell->k);
#else
	if ((cell->k - 1) % 100 == 0) {
		printf("\nIteration-%4d: ", cell->k);
	}
#endif

	cell->candidEst = vXvSparse(cell->candidX, prob[0]->dBar) + maxCutHeight(cell->cuts, cell->sampleSize, cell->candidX, prob[0]->num->cols, prob[0]->lb);

#if defined(LPMIP_PRINT)
	printf("\ninside callback before SD\n");
	printVector(cell->candidX, cell->master->mac - 1, NULL);
	printf("\nEstimate %0.4f  -  Candidate Est %0.4f\n\n", cell->incumbEst, cell->candidEst);
#endif // defined(LPMIP_PRINT)

	if (isCandidInt(cell->candidX, prob[0]->num->cols))
	{

		if (IPoptimal(prob, cell))
		{
			(*breakLoop) = true; return 0;
		}

		/******* 1. Generate new observations, and add it to the set of observations *******/
		cell->sampleSize += config.SAMPLE_INCREMENT;
		for (obs = 0; obs < config.SAMPLE_INCREMENT; obs++) {
			/* (a) Use the stoc file to generate observations */
			generateOmega(stoc, observ, config.TOLERANCE, &config.RUN_SEED[0], NULL);

			/* (b) Since the problem already has the mean values on the right-hand side, remove it from the original observation */
			for (m = 0; m < stoc->numOmega; m++)
				observ[m] -= stoc->mean[m];

			/* (d) update omegaType with the latest observation. If solving with incumbent then this update has already been processed. */
			cell->sample->omegaIdx[obs] = calcOmega(observ - 1, 0, prob[1]->num->numRV, cell->omega, &cell->sample->newOmegaFlag[obs], config.TOLERANCE);
		}

		/******* 2. Solve the subproblem with candidate solution, form and update the candidate cut *******/
		if ((candidCut = formSDCut(prob, cell, cell->candidX, prob[0]->lb, 1)) < 0) {
			errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
			return 1;
		}

		/******* 3. Solve subproblem with incumbent solution, and form an incumbent cut *******/
		//REmove condition
		if (((cell->k - cell->iCutUpdt) % config.TAU == 0)) {
			if ((cell->iCutIdx = formSDCut(prob, cell, cell->incumbX, prob[0]->lb, 1)) < 0) {
				errMsg("algorithm", "solveCell", "failed to create the incumbent cut", 0);
				return 1;
			}
			cell->iCutUpdt = cell->k;
		}

		/******* 4. Check improvement in predicted values at candidate solution *******/
		if (true)
			/* If the incumbent has not changed in the current iteration */
			checkImprovement_callback(prob[0], cell, candidCut);

		/******* 5. Solve the master problem to obtain the new candidate solution */
		if (solveLPMaster(prob[0]->num, prob[0]->dBar, cell, prob[0]->lb)) {
			errMsg("algorithm", "solveCell", "failed to solve master problem", 0);
			return 1;
		}


	}
	else
	{


		/******* 1. Generate new observations, and add it to the set of observations *******/
		cell->sampleSize += config.SAMPLE_INCREMENT;
		for (obs = 0; obs < config.SAMPLE_INCREMENT; obs++) {
			/* (a) Use the stoc file to generate observations */
			generateOmega(stoc, observ, config.TOLERANCE, &config.RUN_SEED[0], NULL);

			/* (b) Since the problem already has the mean values on the right-hand side, remove it from the original observation */
			for (m = 0; m < stoc->numOmega; m++)
				observ[m] -= stoc->mean[m];

			/* (d) update omegaType with the latest observation. If solving with incumbent then this update has already been processed. */
			cell->sample->omegaIdx[obs] = calcOmega(observ - 1, 0, prob[1]->num->numRV, cell->omega, &cell->sample->newOmegaFlag[obs], config.TOLERANCE);
		}

		/******* 2. Solve the subproblem with candidate solution, form and update the candidate cut *******/
		if ((candidCut = formSDCut(prob, cell, cell->candidX, prob[0]->lb, 1)) < 0) {
			errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
			return 1;
		}

		/******* 3. Solve subproblem with incumbent solution, and form an incumbent cut *******/
		//REmove condition
		if (((cell->k - cell->iCutUpdt) % config.TAU == 0)) {
			if ((cell->iCutIdx = formSDCut(prob, cell, cell->incumbX, prob[0]->lb, 1)) < 0) {
				errMsg("algorithm", "solveCell", "failed to create the incumbent cut", 0);
				return 1;
			}
			cell->iCutUpdt = cell->k;
		}

		/******* 4. Check improvement in predicted values at candidate solution *******/
		if (true)
			/* If the incumbent has not changed in the current iteration */
			checkImprovement_callback(prob[0], cell, candidCut);

		/******* 5. Solve the master problem to obtain the new candidate solution */
		if (solveLPMaster(prob[0]->num, prob[0]->dBar, cell, prob[0]->lb)) {
			errMsg("algorithm", "solveCell", "failed to solve master problem", 0);
			return 1;
		}

		///******* 6. Optimality tests *******/
		if (LPoptimal(prob, cell))
		{
			(*breakLoop) = true; return 0;
		}

	}


#if defined(LPMIP_PRINT)
	printf("\ninside callback after SD \n");
	printVector(cell->incumbX, cell->master->mac - 1, NULL);
	printf("\nEstimate %0.4f\n\n", cell->incumbEst);
#endif // defined(LPMIP_PRINT)


	return 0;

}


/*
  Turn the QP master problem to LP
  Siavash Tabrizian May 20
*/
int QPtoLP(stocType *stoc, probType **prob, cellType *cell, int toMIP) {

	cString lu; cString uu; 
	iVector indices;
	int numCols = prob[0]->num->cols;
	int status = 0;

	if (!(indices = (iVector)arr_alloc(numCols, int)))
		errMsg("allocation", "bendersCallback", "indices", 0);
	for (int c = 0; c < numCols; c++)
		indices[c] = c;

	if (toMIP == 0)
	{
	
		/*********00. Set sigma to 0 *********/
		if (changeQPproximal(cell->master->lp, prob[0]->num->cols, 0.0)) {
			errMsg("algorithm", "SolveIntCell", "failed to change the proximal term", 0);
			return 1;
		}

		/*********01. Update the right hand sides of master *********/
		if (revchangeQPrhs(prob[0], cell, cell->incumbX)) {
			errMsg("algorithm", "SolveIntCell", "failed to change the right-hand side to convert the problem into QP", 0);
			return 1;
		}

		/*********02. Change LP solver to primal simplex *********/
		if (changeProbType(cell->master->lp, PROB_LP))
		{
			errMsg("algorithm", "SolveIntCell", "failed to set primal algorithm for LP", 0);
			return 1;
		}
		if (changeLPSolverType(PROB_LP))
		{
			errMsg("algorithm", "SolveIntCell", "failed to set primal algorithm for LP", 0);
			return 1;
		}

	}
	else
	{
		/*********00. Set sigma to 0 *********/
		if (changeQPproximal(cell->master->lp, prob[0]->num->cols, 0.0)) {
			errMsg("algorithm", "SolveIntCell", "failed to change the proximal term", 0);
			return 1;
		}

		/*********01. Update the right hand sides of master *********/
		if (revchangeQPrhs(prob[0], cell, cell->incumbX)) {
			errMsg("algorithm", "SolveIntCell", "failed to change the right-hand side to convert the problem into QP", 0);
			return 1;
		}

		cell->master->type = PROB_MILP;
		
		if (changeCtype(cell->master->lp, numCols, indices, cell->master->ctype)) {
			errMsg("solver", "bendersCallback", "failed to change column type", 0);
			return 1;
		}

		if (changeProbType(cell->master->lp, PROB_MILP)) {
			errMsg("Problem Setup", "bendersCallback", "master", 0);
			return 0;
		}

		/*********02. Change LP solver to B&B *********/
		if (changeLPSolverType(PROB_MILP))
		{
			errMsg("algorithm", "SolveIntCell", "failed to set primal algorithm for LP", 0);
			return 1;
		}

	}
	prob[0]->sp->bdl;
	///* Change the column types to integer/binary to change the problem type to SMIP */
	if (!(lu = (cString)arr_alloc(numCols, int)))
		errMsg("allocation", "algo-copymaster", "lu", 0);
	if (!(uu = (cString)arr_alloc(numCols, int)))
		errMsg("allocation", "algo-copymaster", "uu", 0);


	for (int c = 0; c < numCols; c++)
	{
		lu[c] = 'L';
		uu[c] = 'U';
	}
	status = changeBDS(cell->master->lp, numCols, indices, lu, prob[0]->sp->bdl);
	if (status) {
		solverErrmsg(status);
		errMsg("solver", "algo-copymaster", "failed to change column lowerbound", 0);
		return 1;
	}
	status = changeBDS(cell->master->lp, numCols, indices, uu, prob[0]->sp->bdu);
	if (status) {
		solverErrmsg(status);
		errMsg("solver", "algo-copymaster", "failed to change column upperbound", 0);
		return 1;
	}


#if defined(LPMIP_PRINT)
	if (toMIP == 0)
	{
		writeProblem(cell->master->lp, "finalMaster_QP2LP.lp");
	}
	else
	{
		writeProblem(cell->master->lp, "finalMaster_QP2MILP.lp");
	}
#endif

	mem_free(indices);


}

/*
Turn the LP master problem to MILP
Siavash Tabrizian July 20
*/
int LPtoMILP(stocType *stoc, probType **prob, cellType *cell) {

	cString lu; cString uu;
	iVector indices;
	int numCols = prob[0]->num->cols;
	int status = 0;

	if (!(indices = (iVector)arr_alloc(numCols, int)))
		errMsg("allocation", "bendersCallback", "indices", 0);
	for (int c = 0; c < numCols; c++)
		indices[c] = c;

	cell->master->type = PROB_MILP;

	if (changeCtype(cell->master->lp, numCols, indices, cell->master->ctype)) {
		errMsg("solver", "bendersCallback", "failed to change column type", 0);
		return 1;
	}

	if (changeProbType(cell->master->lp, PROB_MILP)) {
		errMsg("Problem Setup", "bendersCallback", "master", 0);
		return 0;
	}

	/*********02. Change LP solver to B&B *********/
	if (changeLPSolverType(PROB_MILP))
	{
		errMsg("algorithm", "SolveIntCell", "failed to set primal algorithm for LP", 0);
		return 1;
	}


#if defined(LPMIP_PRINT)
	writeProblem(cell->master->lp, "finalMaster_LP2MILP.lp");
#endif

	mem_free(indices);


}

/*
Copy the cell to a cloned cell
Siavash Tabrizian July 20
*/
int copyCell(cellType *cell, cellType *clone_cell)
{
	*clone_cell = *cell;
	clone_cell->master->type = MILP;
	LPptr masterLP;

	/* Copy the master problem as a SMIP into another LPptr associated with the new environment. */
	if (copyMasterSMIP(getEnv(), &masterLP, cell, cell->master->mac)) {
		errMsg("algorithm", "bendersCallback", "failed to create a copy of the master problem", 0);
		return 1;
	}

	clone_cell->master->lp = masterLP;


	return 0;
}

/* This subroutine copies all the elements of the problem from the master problem in the cell structure. This includes the original
* master problem as well as the cuts that were generated in the first-phase of the algorithm */
int copyMasterSMIP(ENVptr env, LPptr *lp, cellType *cell, int numCols) {
	//iVector indices;
	//cString lu; cString uu; 
	//dVector lb; dVector ub;
	int status;

	/* Create a clone of the current master problem in the cell structure, albeit with the new environment. */
	(*lp) = CPXcloneprob(env, cell->master->lp, &status);
	if (status) {
		solverErrmsg(status);
		return 1;
	}


	///* Change the column types to integer/binary to change the problem type to SMIP */
	//if (!(indices = (iVector)arr_alloc(numCols, int)))
	//	errMsg("allocation", "algo-copymaster", "indices", 0);
	//if (!(lu = (cString)arr_alloc(numCols, int)))
	//	errMsg("allocation", "algo-copymaster", "lu", 0);
	//if (!(uu = (cString)arr_alloc(numCols, int)))
	//	errMsg("allocation", "algo-copymaster", "uu", 0);
	//if (!(lb = (dVector)arr_alloc(numCols, int)))
	//	errMsg("allocation", "algo-copymaster", "lb", 0);
	//if (!(ub = (dVector)arr_alloc(numCols, int)))
	//	errMsg("allocation", "algo-copymaster", "ub", 0);

	//for (int c = 1; c < numCols - 1; c++)
	//{
	//	indices[c] = c;
	//	lu[c] = 'L';
	//	uu[c] = 'U';
	//	lb[c] = 0.0;
	//	ub[c] = 1.0;
	//}
	//status = changeBDS((*lp), numCols - 1, indices, lu, lb);
	//if (status) {
	//	solverErrmsg(status);
	//	errMsg("solver", "algo-copymaster", "failed to change column lowerbound", 0);
	//	return 1;
	//}
	//status = changeBDS((*lp), numCols - 1, indices, uu, ub);
	//if (status) {
	//	solverErrmsg(status);
	//	errMsg("solver", "algo-copymaster", "failed to change column upperbound", 0);
	//	return 1;
	//}
	//status = CPXchgctype(env, (*lp), numCols, indices, cell->master->ctype);
	//if (status) {
	//	solverErrmsg(status);
	//	errMsg("solver", "algo-copymaster", "failed to change column type", 0);
	//	return 1;
	//}
	//mem_free(indices);

	///* Change the problem type on the solver. */
	//status = CPXchgprobtype(env, (*lp), PROB_MILP);
	//if (status) {
	//	solverErrmsg(status);
	//	errMsg("Problem Setup", "bendersCallback", "could not change the type of master problem", 0);
	//	return 1;
	//}


	return 0;
}//END copyMasterSMIP()

/* Check if the sloution is integer*/
int isCandidInt(dVector candidate, int size)
{
	int output = 1;

	for (int s = 0; s < size; s++)
	{
		if (candidate[s] != floor(candidate[s]))
		{
			output = 0;
			break;
		}
	}

	return output;
}

/* Make a rounding solution from a fractional */
dVector roundX(dVector a, int len)
{
	int		i;
	dVector	b;

	if ((b = (dVector)arr_alloc(len, double))) {
		for (i = 0; i < len; i++)
			b[i] = ceil(a[i]);
	}
	else
		errMsg("algo", "rounding solutioin", "b", 1);

	return b;
}

int solveIntCell(stocType *stoc, probType **prob, cellType *cell) {
	int			m, GMICut, MIRCut;
	clock_t		tic;

	tic = clock();
	MIRCut = 0;
	GMICut = 0;
	
	//if (cell->master->type == PROB_QP)
	//{
	//	writeProblem(cell->master->lp, "QPMaster.lp");

	//	QPtoLP(stoc, prob, cell, 0);

	//	/*********03. solve new master LP *********/
	//	int 	status;
	//	if (solveProblem(cell->master->lp, cell->master->name, PROB_LP, cell->master->mar, cell->master->mac, &status, config.TOLERANCE)) {
	//		if (status == STAT_INFEASIBLE) {
	//			errMsg("solveIntCell", "solveMaster", "Master problem is infeasible. Check the problem formulation!", 0);
	//			writeProblem(cell->master->lp, "infeasibleM.lp");
	//		}
	//		else {
	//			writeProblem(cell->master->lp, "errorM.lp");
	//			errMsg("solveIntCell", "solveMaster", "failed to solve the master problem", 0);
	//		}
	//	}

	//	writeProblem(cell->master->lp, "LPMaster.lp");
	//	printf("Problem type after: %i \n", getProbType(cell->master->lp));
	//	/* Get the most recent optimal solution to master program */
	//	if (getPrimal(cell->master->lp, cell->candidX, prob[0]->num->cols)) {
	//		errMsg("solveIntCell", "solveMaster", "failed to obtain the primal solution for master", 0);
	//		return 1;
	//	}
	//}

	//////***********************************************
	///      Adding MIP cuts to the master problem
	//////***********************************************
//	while (cell->MIPFlag == false)
//	{

		/******* 1. Get the basis, and form a GMI incumbent cut *******/
#if defined(GMIcutsActive)
		if ((GMICut = formGMICut(prob, cell, cell->candidX, prob[0]->lb)) < 0) {
			errMsg("algorithm", "solveCell", "failed to create the GMI incumbent cut", 0);
			goto TERMINATE;
		}
#endif // GMIcutsActive



		/******* 2. Form a MIR incumbent cut *******/
#if defined(MIRcutsActive)
		if ((MIRCut = formMIRCut(prob, cell, cell->incumbX, prob[0]->lb)) < 0) {
			errMsg("algorithm", "solveCell", "failed to create the MIR incumbent cut", 0);
			goto TERMINATE;
		}
#endif // defined(MIRcutsActive)


		/******* 3. Solve the master problem to obtain the new candidate solution */
		cell->gk += GMICut;
		cell->mk += MIRCut;
		
		int 	status;
		if (solveProblem(cell->master->lp, cell->master->name, PROB_LP, cell->master->mar, cell->master->mac, &status, config.INT_TOLERANCE)) {
			if (status == STAT_INFEASIBLE) {
				errMsg("solveIntCell", "solveMaster", "Master problem is infeasible. Check the problem formulation!", 0);
				writeProblem(cell->master->lp, "infeasibleM.lp");
			}
			else {
				writeProblem(cell->master->lp, "errorM.lp");
				errMsg("solveIntCell", "solveMaster", "failed to solve the master problem", 0);
			}
		}

		/* Get the most recent optimal solution to master program */
		if (getPrimal(cell->master->lp, cell->candidX, prob[0]->num->cols)) {
			errMsg("solveIntCell", "solveMaster", "failed to obtain the primal solution for master", 0);
			return 1;
		}


		bool IPflag = true;
		for (int n = 0; n < prob[0]->num->cols; n++)
		{
			if (cell->candidX[n] - floor(cell->candidX[n]) <= config.INT_TOLERANCE)
			{
				IPflag = false;
			}
		}
		cell->MIPFlag == IPflag;

		cell->time.masterAccumTime += cell->time.masterIter; cell->time.subprobAccumTime += cell->time.subprobIter;
		cell->time.argmaxAccumTime += cell->time.argmaxIter; cell->time.optTestAccumTime += cell->time.optTestIter;
		cell->time.masterIter = cell->time.subprobIter = cell->time.optTestIter = cell->time.argmaxIter = 0.0;
		cell->time.iterTime = ((double)clock() - tic) / CLOCKS_PER_SEC; cell->time.iterAccumTime += cell->time.iterTime;

//	}

	return 0;

TERMINATE:
	return 1;
}//END solveCell()

void printNodeInfo(nodeInfo    *nodeSol, int Nodecnt)
{
	printf("\n --------------- Printing the node solutions Info: \n");
	printf("\n number of nodes: %i \n", Nodecnt);
	for (int i = 0; i < Nodecnt; i++)
	{
		printf("\n ----- ------- ----- \n");
		printf("\n call: %i \n", nodeSol[i].nodeNum);
		printVector(nodeSol[i].sol, nodeSol[i].sol_size, NULL);
		printf("\n Estimate: %0.4f", nodeSol[i].LB);
		if (nodeSol[i].isInt)
		{
			printf("\n integer solution!\n");
		}
		printf("\n dual of the master for this node: \n");
		printVector(nodeSol[i].piM, nodeSol[i].mar, NULL);
		printf("\n ----- ------- ----- \n");
	}
	printf("\n ---------------------------------------------------\n");
}

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



