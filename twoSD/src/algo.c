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

#include "twoSD.h"

extern cString outputDir;
extern configType config;


int algo(oneProblem *orig, timeType *tim, stocType *stoc, cString inputDir, cString probName) {
	dVector	 meanSol = NULL;
	probType **prob = NULL;
	cellType *cell = NULL;
	dVector	lb = NULL;
	batchSummary *batch = NULL;
	FILE 	*sFile = NULL, *iFile = NULL, *dFile = NULL, *oFile = NULL;

	/* open solver environment */
	openSolver();

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell, &batch, &meanSol, lb) )
		goto TERMINATE;

	cell->meanVal = orig->objective;

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

#if  defined(experiment)
		/* evaluating a test solution*/
		double eval2;
		double candidX[11] = { 30.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0 };
		evaluate(sFile, stoc, prob, cell->subprob, candidX);
		goto TERMINATE;
#endif //  defined(experiment)

		/* Use two-stage stochastic decomposition algorithm to solve the problem */
		if (branchbound(stoc, prob, cell, -INFINITY, INFINITY)) {
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
		writeOptimizationStatistics(dFile, iFile, prob, cell, rep);

		/* evaluate the optimal solution*/
		if (config.EVAL_FLAG == 1)
		{
			//LB estimate: CTx (dot product) + height of the cut(subroutine: )
			//cell->incumbEst = vXvSparse(cell->incumbX, prob[0]->dBar) + maxCutHeight(cell->cuts, cell->sampleSize, cell->incumbX, prob[0]->num->cols, cell->lb);
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
	mem_free(lb);
	return 1;
}//END algo()

int solveCell(stocType *stoc, probType **prob, cellType *cell) {
	dVector 	observ;
	clock_t		tic;
	bool breakLoop = false;



	/* -+-+-+-+-+-+-+-+-+-+-+-+-+-+- Main Algorithm -+-+-+-+-+-+-+-+-+-+-+-+-+-+- */
	observ = (dVector) arr_alloc(stoc->numOmega + 1, double);

	/******* 0. Initialization: The algorithm begins by solving the master problem as a QP *******/
	while (cell->optFlag == false && cell->ki < config.MAX_ITER_CLBK && cell->k < config.MAX_ITER) {
		
		tic = clock();

		if (mainloopSDCell(stoc, prob, cell, &breakLoop, observ)) {
			errMsg("Callback", "usersolve", "failed to solve Benders cell for the node problem", 0);
			goto TERMINATE;
		}

		cell->time.masterAccumTime += cell->time.masterIter; cell->time.subprobAccumTime += cell->time.subprobIter;
		cell->time.argmaxAccumTime += cell->time.argmaxIter; cell->time.optTestAccumTime += cell->time.optTestIter;
		cell->time.masterIter = cell->time.subprobIter = cell->time.optTestIter = cell->time.argmaxIter = 0.0;
		cell->time.iterTime = ((double) clock() - tic)/CLOCKS_PER_SEC; cell->time.iterAccumTime += cell->time.iterTime;
		if (cell->ki > 500 && cell->ki % 100 == 0) {
			printf("."); fflush(stdout);
		}

		if (breakLoop)
		{
			break;
		}

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
	cell->ki++;
#if defined(STOCH_CHECK) || defined(ALGO_CHECK)
	printf("\nIteration-%d :: \n", cell->k);
#else
	//if ((cell->k - 1) % 100 == 0) {
	//	printf("\nIteration-%4d: ", cell->k);
	//}
#endif

	if(cell->ki > 1)
	/******* 1. Optimality tests *******/
	if (IPoptimal(prob, cell))
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
	if (cell->ki > 1)
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
	evaluate(NULL,NULL, stoc, prob, cell->subprob, cell->candidX);
	

	//1a - Change the master to MILP no callback solve to optimality - > xIP (LB)
	//1b - check the optimality of (xIP) using pretest 
	//1c - Evaluate the solution (xIP) using evaluate(sFile, stoc, prob, cell->subprob, cell->incumbX); -> (UB)
	/* Turn the clone problem to LP */
	QPtoLP(stoc, prob, cell, 0);

	/* Launch the solver to solve the MIP master problem */
	if (solveProblem(cell->master->lp, cell->master->name, cell->master->type, &status)) {
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
	evaluate(NULL,NULL, stoc, prob, cell->subprob, cell->candidX);

	/* Turn the clone problem to MILP */
	LPtoMILP(stoc, prob, cell);
	/* Launch the solver to solve the MIP master problem */
	if (solveProblem(cell->master->lp, cell->master->name, cell->master->type, &status)) {
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
	evaluate(NULL,NULL, stoc, prob, cell->subprob, cell->candidX);

	//summary for the second one 

	printf("----------------------------------------------------------------------------------------------------\n");

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
		//if (changeQPproximal(cell->master->lp, prob[0]->num->cols, 0.0)) {
		//	errMsg("algorithm", "SolveIntCell", "failed to change the proximal term", 0);
		//	return 1;
		//}

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
			return 1;
		}

		/*********02. Change LP solver to B&B *********/
		if (changeLPSolverType(PROB_MILP))
		{
			errMsg("algorithm", "SolveIntCell", "failed to set primal algorithm for LP", 0);
			return 1;
		}

	}

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

	return 0;
}

/*
Turn the LP master problem to MILP
Siavash Tabrizian July 20
*/
int LPtoMILP(stocType *stoc, probType **prob, cellType *cell) {
	iVector indices;
	int numCols = prob[0]->num->cols;

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
		return 1;
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
	return 0;
}//END LPtoMILP();

/*
Getting the row names in the bnb after cuts are added 
Siavash Tabrizian Sep 20
*/
void getRowNameMaster(cellType *cell)
{
	int numRow = getNumRows(cell->master->lp);
	int cur_rownamespace;
	char          **cur_rowname = NULL;
	char          *cur_rownamestore = NULL;
	int           surplus;
	int			  status;
	cell->rownum = numRow;

	status = getRowName(cell->master->lp, 0, numRow, NULL, NULL, 0, &surplus);

	if ((status != CPXERR_NEGATIVE_SURPLUS) &&
		(status != 0)) {
		fprintf(stderr,
			"Could not determine amount of space for row names.\n");
	}

	cur_rownamespace = -surplus;
	if (cur_rownamespace > 0) {
		cur_rowname = (char **)malloc(sizeof(char *)*numRow);
		cur_rownamestore = (char *)malloc(cur_rownamespace);
		if (cur_rowname == NULL ||
			cur_rownamestore == NULL) {
			fprintf(stderr, "Failed to get memory for row names.\n");
			status = -1;
		}
		status = getRowName(cell->master->lp, 0, numRow, cur_rowname, cur_rownamestore, cur_rownamespace, &surplus);
		if (status) {
			fprintf(stderr, "CPXgetcolname failed.\n");
		}
	}
	else {
		printf("No names associated with problem.  Using Fake names.\n");
	}

	cell->cur_rowname = cur_rowname;
}

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



