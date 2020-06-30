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

int algo(oneProblem *orig, timeType *tim, stocType *stoc, cString inputDir, cString probName) {
	dVector	 meanSol = NULL;
	probType **prob = NULL;
	cellType *cell = NULL;
	batchSummary *batch = NULL;
	FILE 	*sFile = NULL, *iFile = NULL;

	/* open solver environment */
	openSolver();

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell, &batch, &meanSol) )
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
		if ( solveCell(stoc, prob, cell) ) {
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
			printVector(cell->incumbX, cell->master->mac,NULL);
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
	freeCellType(cell);
	freeProbType(prob, 2);
	return 1;
}//END algo()

int intalgo(oneProblem *orig, timeType *tim, stocType *stoc, cString inputDir, cString probName) {
	dVector	 meanSol = NULL;
	probType **prob = NULL;
	cellType *cell = NULL;
	batchSummary *batch = NULL;
	FILE 	*sFile = NULL, *iFile = NULL;

	/* open solver environment */
	openSolver();

	/* complete necessary initialization for the algorithm */
	if (setupAlgo(orig, stoc, tim, &prob, &cell, &batch, &meanSol))
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
		if (solveCell(stoc, prob, cell)) {
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
	return 1;
}//END algo()

int solveCell(stocType *stoc, probType **prob, cellType *cell) {
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

	if (config.SMIP != MILP) {
		/* Phase-1 has completed, we have an approximation obtained by solving the relaxed problem.
		* Phase-2 be used to impose integrality through costom procedures or the callback.  */
		cell->optFlag = false;
		if (bendersCallback(stoc, prob, cell)) {
			errMsg("algorithm", "solverCell", "failed to run the Benders-callback routine", 0);
			goto TERMINATE;
		}
	}

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
	if ((candidCut = formSDCut(prob, cell, cell->candidX, prob[0]->lb)) < 0) {
		errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
		return 1;
	}

	/******* 4. Solve subproblem with incumbent solution, and form an incumbent cut *******/
	if (((cell->k - cell->iCutUpdt) % config.TAU == 0)) {
		if ((cell->iCutIdx = formSDCut(prob, cell, cell->incumbX, prob[0]->lb)) < 0) {
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
	if ((candidCut = formSDCut(prob, cell, cell->candidX, prob[0]->lb)) < 0) {
		errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
		return 1;
	}

	/******* 4. Solve subproblem with incumbent solution, and form an incumbent cut *******/
	if (((cell->k - cell->iCutUpdt) % config.TAU == 0)) {
		if ((cell->iCutIdx = formSDCut(prob, cell, cell->incumbX, prob[0]->lb)) < 0) {
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
	if (solveLPMaster(prob[0]->num, prob[0]->dBar, cell, prob[0]->lb)) {
		errMsg("algorithm", "solveCell", "failed to solve master problem", 0);
		return 1;
	}

	return 0;

}

int QPtoLP(stocType *stoc, probType **prob, cellType *cell, int toMIP) {

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
		iVector indices;

		if (!(indices = (iVector)arr_alloc(prob[0]->num->cols, int)))
			errMsg("allocation", "bendersCallback", "indices", 0);
		for (int c = 0; c < prob[0]->num->cols; c++)
			indices[c] = c;
		if (changeCtype(cell->master->lp, prob[0]->num->cols, indices, cell->master->ctype)) {
			errMsg("solver", "bendersCallback", "failed to change column type", 0);
			return 1;
		}
		mem_free(indices);

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



