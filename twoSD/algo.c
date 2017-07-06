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

extern configType config;
extern int numAgents;
extern string outputDir;

int algo(oneProblem *orig, timeType *tim, stocType *stoc, string inputDir, string probName){
	vector          xk = NULL, lb = NULL;
	probType        **prob = NULL;
	cellType        *cell = NULL;
	clock_t         tic, toc;
	double          totRunTime;

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell) )
		goto TERMINATE;

	/* create cells used in the algorithm */
	cell = newCell(prob, xk, weight);
	if (cell == NULL) {
		errMsg("setup", "algo", "failed to setup the cell used by the algorithm", 0);
		goto TERMINATE;
	}

	tic = clock();
	/* Use MASP algorithm to solve the problem */
	if ( solveMASP(stoc, prob, cell, sim, inputDir, probName) ) {
		errMsg("algorithm", "algo", "failed to solve the cells using MASP algorithm", 0);
		goto TERMINATE;
	}

	toc = clock();
	totRunTime = ((double)(toc - tic))/CLOCKS_PER_SEC;

	writeStat(prob[0], cell, totRunTime, probName);

	if(xk) mem_free(xk);
	if(lb) mem_free(lb);
	if(cell) freeCellType(cell);
	if(prob) freeProbType(prob, numAgents);
	return 0;

	TERMINATE:
	if(xk) mem_free(xk);
	if(lb) mem_free(lb);
	if(cell) freeCellType(cell);
	if(prob) freeProbType(prob, numAgents);
	return 1;
}//END algo()

int solveMASP(stocType *stoc, probType **prob, cellType **cell, simType *sim, string inputDir, string probName) {
	vector 	observ;
	int		status = 0, i, m, sOffset;
	clock_t tic, toc, start, end;
	FILE    *StatTime;

	StatTime = openFile(outputDir, "time.dat", "w");

	/* -+-+-+-+-+-+-+-+-+-+-+-+-+-+- Main Algorithm -+-+-+-+-+-+-+-+-+-+-+-+-+-+- */
	if ( !(observ = (vector) arr_alloc(stoc->numOmega + 1, double)) )
		errMsg("allocation", "solveMASP", "observ", 0);

	/******* 0. Initialization: The algorithm begins by solving the master problem as a QP *******/
	while (cell[0]->optFlag == FALSE && cell[0]->k < config.MAX_ITERATION) {
		start = clock();

		cell[0]->k++;

#if 0
		printf("\nIteration-%d :: \n", cell[0]->k);
#else
		if ( (cell[0]->k -1)% 50 == 0)
			printf("\nIteration-%4d: ", cell[0]->k);
#endif

		/******* 1. Optimality tests *******/
		if (optimal(prob, cell))
			break;

		/******* 2. Generate new observation *******/
		if (config.STOCFILE == 1)
			/* use the stoc file to generate observations */
			generateOmega(stoc, observ, &config.RUN_SEED);
		else
			/* use the simulator to generate observations */
			simulateOmega(sim, observ);
		/* Since the problem already has the mean values on the right-hand side, remove it from the original observation */
		for ( m = 0; m < stoc->numOmega; m++ )
			observ[m] -= stoc->mean[m];

		/******* 3. solve the subproblem with current observation and candidate solution first, and then with incumbent solution *******/
		sOffset = 0;
		for (i = 1; i < numAgents; i++) {
			tic = clock();

			status = solveAgents(cell[i], prob[i], cell[0]->candidX, cell[0]->incumbX, prob[0]->num->cols, observ + sOffset, prob[0]->lb, cell[0]->k);

			if ( status ) {
				errMsg("algorithm", "solveMASP", "failed to solve agent problem", 0);
				return 1;
			}
			if (!(prob[i]->omegas == NULL))
				sOffset += prob[i]->omegas->numRV;

			toc = clock();
			cell[i]->runTime->iterTime = ((double) (toc-tic)) / CLOCKS_PER_SEC;
		}

		/******** 5. Add the new cut(s) to master problem ********/
		status = addCuts2Master(cell, prob[0]->num->cols, prob[0]->lb);
		if ( status ) {
			errMsg("algorithm", "solveMASP", "failed to add cuts to the master problem", 0);
			return 1;
		}

		/******* 4. Check improvement in predicted values at candidate solution *******/
		if ( !(cell[0]->incumbChg) && cell[0]->k > 1)
			/* If the incumbent has not changed in the current iteration */
			checkImprovement(prob[0], cell);

		/******* 5. Solve QP master problem *******/
		status = constructQP(cell[0]->sp->lp, prob[0]->num->cols, cell[0]->quadScalar);
		if ( status ) {
			errMsg("algorithm", "solveMASP", "failed to change the proximal term", 0);
			return 1;
		}

		status = solveQPMaster(prob[0]->num, prob[0]->dBar, cell, prob[0]->sp->mar, prob[0]->lb);
		if ( status ) {
			errMsg("algorithm", "solveMASP", "failed to solve master problem", 0);
			return 1;
		}

		end = clock();
		cell[0]->runTime->iterTime = ((double) (end-start)) / CLOCKS_PER_SEC;

		for (i = 0; i < numAgents; i++) {
			cell[i]->runTime->totSolTime += cell[i]->runTime->iterSolTime;
			cell[i]->runTime->totCutGenTime += cell[i]->runTime->iterCutGenTime;
			cell[i]->runTime->totTime += cell[i]->runTime->iterTime;
#ifdef RUNTIME
			if ( i == 0) {
				printf("%d\t", cell[0]->k);
				fprintf(StatTime, "%d\t", cell[0]->k);
			}
			printf("%lf%lf\t%lf", cell[i]->runTime->iterTime, cell[i]->runTime->iterSolTime, cell[i]->runTime->iterCutGenTime);
			fprintf(StatTime, "%lf\t%lf\t%lf\t", cell[i]->runTime->iterTime, cell[i]->runTime->iterSolTime, cell[i]->runTime->iterCutGenTime);
			if ( i == numAgents-1 ) {
				printf("\n");
				fprintf(StatTime, "\n");
			}
#endif
		}
	}//END while loop

	writeStat(prob[0], cell, 0.0, probName);

	/*evaluating the optimal solution*/
	if (config.EVAL_FLAG == 1) {
		if (config.STOCFILE == 1) {
			evaluateOptStoc(stoc, prob, cell);
		}
		else {
			evaluateOptSim(stoc, prob, cell, sim);
		}
	}

	mem_free(observ);
	fclose(StatTime);

	return 0;
}//END solveMASP()

int solveAgents(cellType *cell, probType *prob, vector candidX, vector incumbX, int lenX, vector observ, double lb, int iter) {
	int     status,omegaIdx;
	BOOL	newOmegaFlag = FALSE; /* when creating incumbent cut, newOmegaFlag should be false */
	clock_t     tic, toc;

	/* (a) update omegaType with the latest observation. If solving with incumbent then this update has already been processed. */
	omegaIdx = calcOmega(observ - 1, 0, prob->num->numRV, cell->omega, &newOmegaFlag);

	/* (b) obtain recourse information form the agents. If the agent is turned ON then a subproblem is solved and the new dual information is obtained.
	 * On the other hand, if the subproblem is turned OFF, then previous dual information is presumed to be sufficient, and only updates are made to
	 * the delta structure */
	if ( !(cell->optFlag) ) {
		/* Agent is ON */
		cell->k++;

		status = solveSubprob(prob, cell, candidX, cell->omega->vals[omegaIdx]);
	}

	/* (c) update the stochastic elements in the problem */
	stochasticUpdates(prob->num, candidX, prob->coord, prob->bBar, prob->Cbar, cell->lambda, cell->sigma, cell->delta, cell->omega,
			newOmegaFlag, omegaIdx, config.MAX_ITERATION, cell->k, cell->pi, cell->mubBar, cell->optFlag);

	/* change the new omega flag to reflect that stochastic updates have been completed */
	newOmegaFlag = FALSE;

	/* (d) Generate candidate cut */
	tic = clock();

	cell->cCutIdx = formSDcut(prob, cell, candidX, lenX, lb, iter);

	toc = clock();
	cell->runTime->iterCutGenTime = ((double) (toc-tic)) / CLOCKS_PER_SEC;
	cell->runTime->totCutGenTime += cell->runTime->iterCutGenTime;

	if ( cell->cCutIdx < 0 ) {
		errMsg("algorithm", "solveAgents", "failed to create candidate cut for agent-i", 0);
		return 1;
	}

	/* (e) Solve the problem again with incumbent solution once every TAU iteration */
	if (((iter - cell->iCutUpdt) % config.TAU == 0 ) ) {
		if ( !(cell->optFlag) ) {
			status = solveSubprob(prob, cell, incumbX, cell->omega->vals[omegaIdx]);
			if ( status ) {
				errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
				return 1;
			}
		}

		/* (f) Stochastic updates with respect to incumbent information */
		stochasticUpdates(prob->num, candidX, prob->coord, prob->bBar, prob->Cbar, cell->lambda, cell->sigma, cell->delta, cell->omega,
				newOmegaFlag, omegaIdx, config.MAX_ITERATION, cell->k, cell->pi, cell->mubBar, cell->optFlag);

		/* (f) Generate incumbent cut */
		tic = clock();

		if ( config.MULTI_CUT)
			status = formSDcut(prob, cell, incumbX, lenX, lb, iter);
		else{
			cell->iCutIdx = formSDcut(prob, cell, incumbX, lenX, lb, iter);
			status = cell->iCutIdx;
		}

		toc = clock();

		cell->runTime->iterCutGenTime += ((double) (toc-tic)) / CLOCKS_PER_SEC;

		if ( status < 0 ) {
			errMsg("algorithm", "solveAgents", "failed to create incumbent cut for agent-i", 0);
			return 1;
		}

		cell->iCutUpdt = iter;
	}

	return 0;
}//END solveAgents

/* This function is used to create cells used in the algorithm */
cellType **newCell(probType **prob, vector xk, vector weight) {
	cellType        **cell;
	double   AggWeight;
	int agentCnt;

	/* allocate memory to all cells used in the algorithm. The first cell belongs to the master problem, while the rest correspond to each of the
	 * sub-agents in the problem.  */
	if (!(cell = (cellType **) arr_alloc (numAgents, cellType *)))
		errMsg("Memory allocation", "new_cell", "failed to allocate memory to cell",0);

	cell->master = newMaster(prob[0], xk, weight, AggWeight);
	if ( cell->master == NULL ) {
		errMsg("setup", "newCell", "failed to setup the master problem", 0);
		return NULL;
	}



	/* setup subproblem cells */
	for (agentCnt = 1; agentCnt < numAgents; agentCnt++) {
		AggWeight += weight[agentCnt];
		cell[agentCnt] = newSubprob(prob[agentCnt], agentCnt, weight[agentCnt]);
		if ( cell[agentCnt] == NULL ) {
			errMsg("setup", "newCell", "failed to setup the subproblem cell", 0);
			return NULL;
		}
	}



	return cell;
}//END newCell()

void writeStat(probType *prob, cellType **cell, double totRunTime, string probName) {
	double  Est = 0.0;
	int     i, temp;
	FILE    *Stat;
	FILE    *StatTime;

	StatTime = openFile(outputDir, "time.dat", "a");
	for (i = 0; i < numAgents; i++) {
		if ( i == 0)
			fprintf(StatTime, "\nTotal (%d):\t", cell[0]->k-1);
		fprintf(StatTime, "%lf\t%lf\t%lf\t", cell[i]->runTime->totTime, cell[i]->runTime->totSolTime, cell[i]->runTime->totCutGenTime);
	}

	fprintf(StatTime, "\ntotRunTime = %lf", totRunTime);

	fclose(StatTime);

	Stat = openFile(outputDir, "statistic.dat", "w");
	fprintf(Stat, "-------- Result summary for %s --------\n", probName);

	if ( config.STOCFILE )
		fprintf(Stat, "Samples simulated using SMPS stoch file.\n");
	else
		fprintf(Stat, "Using samples simulated outside the code.\n");

	if ( config.MULTI_CUT )
		fprintf(Stat, "Multiple cuts, each corresponding to individual agents, added to master problem.\n");
	else
		fprintf(Stat, "A single aggregated cut added to master problem.\n\n");

	for (i = 1; i < numAgents; i++) {
		fprintf(Stat, "Agent %d ::\n\tNumber of linear programs = %d\n", cell[i]->ID, cell[i]->k);
		if (config.MULTI_CUT) {
			fprintf(Stat, "\tEstimate at the incumbent solution = %lf\n", cell[i]->incumbEst);
		}
		else{
			Est = maxCutHeight(cell[0]->lbType, cell[i]->cuts, cell[0]->k-1, cell[0]->incumbX, prob->num->cols, prob->lb, &temp);
			fprintf(Stat, "\tEstimate at the incumbent solution = %lf\n", Est);
		}
	}

	fprintf(Stat, "\nTotal number of iteration: %d\n", cell[0]->k-1);
	fprintf(Stat, "Aggregated estimate at the incumbent solution = %lf\n", cell[0]->incumbEst);
	fprintf(Stat, "\nIncumbent solution = ");
	printVector(cell[0]->incumbX, prob->num->cols, Stat);

	fclose(Stat);
}//END WriteStat

void freeCellType(cellType **cell) {
	int n;

	for ( n = 0; n < numAgents; n++ ) {
		if ( cell[n] ) {
			if (n == 0) freeOneProblem(cell[n]->sp);
			if (cell[n]->pi) mem_free(cell[n]->pi);
			if (cell[n]->di) mem_free(cell[n]->di);
			if (cell[n]->cuts) freeCuts(cell[n]->cuts);
			if (cell[n]->fcuts) freeCuts(cell[n]->fcuts);
			if (cell[n]->incumbX) mem_free(cell[n]->incumbX);
			if (cell[n]->candidX) mem_free(cell[n]->candidX);
			if (cell[n]->pi_ratio) mem_free(cell[n]->pi_ratio);
			if (cell[n]->delta) freeDeltaType(cell[n]->delta, cell[n]->lambda->cnt, cell[n]->omega->cnt);
			if (cell[n]->lambda) freeLambdaType(cell[n]->lambda);
			if (cell[n]->sigma) freeSigmaType(cell[n]->sigma);
			if (cell[n]->omega) freeOmegaType(cell[n]->omega);
			if (cell[n]->spRHS) mem_free(cell[n]->spRHS);
			if (cell[n]->runTime) mem_free(cell[n]->runTime);
			mem_free(cell[n]);
		}
	}
	mem_free(cell);

}//END freeCellType()
