/*
 * cell.c
 *
 *  Created on: Feb 24, 2021
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"

extern configType config;

int solveCell(stocType *stoc, probType **prob, cellType *cell) {
	int SDstatus;
	dVector 	observ;
	clock_t		tic;
	int 		candidCut;

	/* -+-+-+-+-+-+-+-+-+-+-+-+-+-+- Main Algorithm -+-+-+-+-+-+-+-+-+-+-+-+-+-+- */
	observ = (dVector) arr_alloc(stoc->numOmega + 1, double);

	/******* 0. Initialization: The algorithm begins by solving the master problem as a QP *******/
	while (cell->optFlag == false && cell->ki < config.MAX_ITER_CLBK && cell->k < config.MAX_ITER) {

		tic = clock();

		/* Main loop of SD */
		cell->k++;
		cell->ki++;
#if defined(STOCH_CHECK) || defined(ALGO_CHECK)
		printf("\nIteration-%d :: \n", cell->k);
#else
		if ((cell->k - 1) % 100 == 0) {
			printf("\nIteration-%4d: ", cell->k);
		}
#endif

		if(cell->ki > 1)
			/******* 1. Optimality tests *******/
			if ( optimal(prob, cell) ) {
				return 0;
			}

		/******* 2. Generate new observations, and add it to the set of observations *******/
		cell->sampleSize += config.SAMPLE_INCREMENT;
		for (int obs = 0; obs < config.SAMPLE_INCREMENT; obs++) {
			/* (a) Use the stoc file to generate observations */
			generateOmega(stoc, observ, config.TOLERANCE, &config.RUN_SEED[0], NULL);

			/* (b) Since the problem already has the mean values on the right-hand side, remove it from the original observation */
			for (int m = 0; m < stoc->numOmega; m++)
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
		if ( solveQPMaster(prob[0]->num, prob[0]->dBar, cell, prob[0]->lb) ) {
			errMsg("algorithm", "solveCell", "failed to solve master problem", 0);
			return 1;
		}

		/* Return when the master problem is infeasible. The infeasibility is handled differently for SLPs and SMIPs. */
		if ( !cell->masterFeasFlag )
			return 1;

		cell->time.masterAccumTime += cell->time.masterIter; cell->time.subprobAccumTime += cell->time.subprobIter;
		cell->time.argmaxAccumTime += cell->time.argmaxIter; cell->time.optTestAccumTime += cell->time.optTestIter;
		cell->time.masterIter = cell->time.subprobIter = cell->time.optTestIter = cell->time.argmaxIter = 0.0;
		cell->time.iterTime = ((double) clock() - tic)/CLOCKS_PER_SEC; cell->time.iterAccumTime += cell->time.iterTime;
		if (cell->ki > 500 && cell->ki % 100 == 0) {
			printf("."); fflush(stdout);
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
	return SDstatus;
}//END solveCell()

int mainloopSDCell(stocType *stoc, probType **prob, cellType *cell, bool *breakLoop, dVector observ) {
	int			m, candidCut, obs;

	cell->k++;
	cell->ki++;
#if defined(STOCH_CHECK) || defined(ALGO_CHECK)
	printf("\nIteration-%d :: \n", cell->k);
#else
	if ((cell->k - 1) % 100 == 0) {
		printf("\nIteration-%4d: ", cell->k);
	}
#endif

	if(cell->ki > 1)
		/******* 1. Optimality tests *******/
		if (IPoptimal(prob, cell)) {
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
	if ( solveQPMaster(prob[0]->num, prob[0]->dBar, cell, prob[0]->lb) ) {
		errMsg("algorithm", "solveCell", "failed to solve master problem", 0);
		return 1;
	}

	if ( !cell->masterFeasFlag )
		(*breakLoop) = true;

	return 0;
}//mainloopSDCell()
