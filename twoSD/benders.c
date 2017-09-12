/*
 * algo.c
 *
 *  Created on: Sep 30, 2015
 *      Author: Harsha Gangammanavar
 */

#include "benders.h"

extern configType config;
extern string outputDir;

int benders (oneProblem *orig, timeType *tim, stocType *stoc, string probName) {
	probType **prob = NULL;
	cellType *cell = NULL;
	FILE 	*soln;

	/* read algorithm configuration file */
	if ( readConfig() )
		goto TERMINATE;

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell) )
		goto TERMINATE;

	/* Update omega structure */
	updateOmega(stoc, cell->omega);

	printf("Starting Benders decomposition.\n");
	/* Use two-stage algorithm to solve the problem */
	if ( solveBendersCell(stoc, prob, cell) ) {
		errMsg("algorithm", "algo", "failed to solve the cells using MASP algorithm", 0);
		goto TERMINATE;
	}

	/* Write solution statistics for optimization process */
	soln = openFile(outputDir, "summary.dat", "w");
	writeBendersStatistic(soln, prob, cell, probName, tim->numStages);
	writeBendersStatistic(stdout, prob, cell, probName, tim->numStages);

	/* evaluating the optimal solution*/
	if (config.EVAL_FLAG == 1)
		evaluateBenders(&soln, stoc, prob, cell, cell->incumbX);

	printf("\nSuccessfully completed the L-shaped method.\n");

	/* free up memory before leaving */
	freeCellType(cell);
	freeProbType(prob, 2);
	return 0;

	TERMINATE:
	if(cell) freeCellType(cell);
	if(prob) freeProbType(prob, 2);
	return 1;
}//END benders()

int solveBendersCell(stocType *stoc, probType **prob, cellType *cell) {
	int 	candidCut;

	/* Main loop of the algorithm */
	while (TRUE) {
		cell->k++;

#if defined(STOCH_CHECK) || defined(ALGO_CHECK)
		printf("\nIteration-%d :: \n", cell->k);
#else
		if ( (cell->k-1) % 100 == 0)
			printf("\nIteration-%4d: ", cell->k);
#endif

		/******* 1. Solve the subproblem with candidate solution, form and update the candidate cut *******/
		if ( (candidCut = formBendersCut(prob[1], cell, cell->candidX, FALSE)) < 0 ) {
			errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
			return 1;
		}
		cell->ub = vXvSparse(cell->candidX, prob[0]->dBar) + cutHeight(cell->cuts->vals[candidCut], cell->candidX, prob[0]->num->cols, FALSE, 0, 0.0);

		/******* 2. Optimality tests *******/
		/* Begin by computing the upper bound, it is the candidate cut height at the candidate solution */
		if (optimalBenders(prob, cell))
			break;

		/******* 3. Solve the master problem to obtain the new candidate solution */
		if ( solveBendersMaster(prob[0]->num, prob[0]->dBar, cell) ) {
			errMsg("algorithm", "solveMASP", "failed to solve master problem", 0);
			return 1;
		}
	}

	return 0;
}//END solveCell()

void updateOmega(stocType *stoc, omegaType *omega) {
	int cnt, i, base, idx;
	BOOL createSAA = FALSE;

	mem_free(omega->weights); omega->weights = NULL;

	if ( strstr(stoc->type, "BLOCKS") != NULL ) {
		if ( (omega->cnt = stoc->numVals[0]) <= config.MAX_OBS) {
			omega->vals = (vector *) mem_realloc(omega->vals, omega->cnt*sizeof(vector));
			if ( !(omega->probs = (vector) arr_alloc(omega->cnt, double)))
				errMsg("allocation", "newOmega", "omega->probs", 0);
			for ( cnt = 0; cnt < omega->cnt; cnt++) {
				omega->probs[cnt]= stoc->probs[0][cnt];
				if ( !(omega->vals[cnt] = (vector) arr_alloc(omega->numRV+1, double)) )
					errMsg("allocation", "newOmega", "omega->vals[cnt]", 0);
				for (i = 0; i < omega->numRV; i++)
					omega->vals[cnt][i+1]=stoc->vals[i][cnt]-stoc->mean[i];
				omega->vals[cnt][0] = oneNorm(omega->vals[cnt]+1, omega->numRV);
			}
		}
		else
			createSAA = TRUE;
	}
	else if ( strstr(stoc->type, "INDEP") != NULL ) {
		omega->cnt = 1;
		for ( i = 0; i < stoc->numOmega; i++ )
			omega->cnt *= stoc->numVals[i];

		if ( omega->cnt > config.MAX_OBS )
			createSAA = TRUE;
		else {
			omega->vals = (vector *) mem_realloc(omega->vals, omega->cnt*sizeof(vector));
			if ( !(omega->probs = (vector) arr_alloc(omega->cnt, double)))
				errMsg("allocation", "newOmega", "omega->probs", 0);
			for ( cnt = 0; cnt < omega->cnt; cnt++) {
				if ( !(omega->vals[cnt] = (vector) arr_alloc(omega->numRV+1, double)) )
					errMsg("allocation", "newOmega", "omega->vals[cnt]", 0);
				omega->probs[cnt] = 1; base = omega->cnt;
				for ( i = 0; i < omega->numRV; i++ ) {
					base /= stoc->numVals[i];
					idx = (int)((double) cnt / (double) base) % stoc->numVals[i];
					omega->vals[cnt][i+1] = stoc->vals[i][idx]-stoc->mean[i];
					omega->probs[cnt] *= stoc->probs[i][idx];
				}
			}
		}
	}
	else {
		printf("\nStoch file has a continuous distribution. Setting up a Sample average approximation.\n");
		createSAA = TRUE;
	}

	if ( config.SAA == 1 && createSAA ) {
		printf("Include procedure to create SAA.\n");
		//			omega->vals = setupSAA(stoc, &config.RUN_SEED, &omega->cnt);
	}

}//END updateOmega()

/* In this his function the master problem is solved after the newest cut is added to master problem, the incumbent cut is updated if necessary.
 * Here the coefficients on all the cuts are updated, and finally master problem is solved. */
int solveBendersMaster(numType *num, sparseVector *dBar, cellType *cell) {
	double 	d2;
	int 	status, i;

#if defined(ALGO_CHECK)
	writeProblem(cell->master->lp,"cellMaster.lp");
#endif

	/* solve the master problem */
	if ( solveProblem(cell->master->lp, cell->master->name, config.MASTER_TYPE, &status) ) {
		writeProblem(cell->master->lp, "error.lp");
		errMsg("algorithm", "solveMaster", "failed to solve the master problem", 0);
		return 1;
	}

	cell->candidEst = getObjective(cell->master->lp, config.MASTER_TYPE);

	/* Get the most recent optimal solution to master program */
	if ( getPrimal(cell->master->lp, cell->candidX, num->cols) ) {
		errMsg("algorithm", "solveMaster", "failed to obtain the primal solution for master", 0);
		return 1;
	}

	/* increment the number of problems solved during algorithm */
	cell->LPcnt++;

	if ( cell->master->type == PROB_QP ) {
		/* Get the dual solution too */
		status = getDual(cell->master->lp, cell->piM, cell->master->mar);
		if ( status ) {
			errMsg("solver", "solveQPMaster", "failed to obtain dual solutions to master", 0);
			return 1;
		}

		/* add the incumbent back to change from \Delta X to X */
		for (i = 1; i <= num->cols; i++)
			d2 += cell->candidX[i] * cell->candidX[i];
		addVectors(cell->candidX, cell->incumbX, NULL, num->cols);

		/* update d_norm_k in soln_type. */
		if (cell->k == 1)
			cell->normDk_1 = d2;
		cell->normDk = d2;

		/* Calculate gamma for next improvement check on incumbent x. */
		cell->gamma = cell->candidEst - cell->incumbEst;
	}

	return 0;
}//END solveMaster()

int formBendersCut(probType *prob, cellType *cell, vector Xvect, BOOL isIncumb) {
	oneCut 	*cut;
	vector 	piCbarX;
	intvec	istar;
	double	multiplier, val, argmax;
	int    	cutIdx, obs, c, cnt, lambdaIdx, sigmaIdx, offset;

	if (!(istar = (intvec) arr_alloc(cell->omega->cnt, int)) )
		errMsg("allocation", "formSDCut", "istar", 0);
	if (!(piCbarX= arr_alloc(cell->sigma->cnt, double)))
			errMsg("Allocation", "SDCut", "pi_Tbar_x",0);
	/* Calculate (Pi x Cbar) x X by mult. each VxT by X, one at a time */
	for (cnt = 0; cnt < cell->sigma->cnt; cnt++)
		piCbarX[cnt] = vXv(cell->sigma->vals[cnt].piC, Xvect, prob->coord->colsC, prob->num->cntCcols);
	offset = prob->num->rvbOmCnt + prob->num->rvCOmCnt;

	/* Only a fraction (at least one) of subproblems are solved in any iteration, for the remainder the argmax operator is used */
	cnt = cell->LPcnt;
	while (cnt == cell->LPcnt) {
		for ( obs = 0; obs < cell->omega->cnt; obs++ ) {
			val = randUniform(&config.SUBPROB_SAMPLE_SEED);

			if ( val <= config.SUBPROB_SAMPLE_PCT ) {
				/* (a) Construct the subproblem with a given observation and master solution, solve the subproblem, and complete stochastic updates. */
				istar[obs] = solveSubprob(prob, cell->subprob, Xvect, cell->basis, cell->lambda, cell->sigma, cell->delta, config.MAX_ITER,
						cell->omega, obs, FALSE, cell->k, config.TOLERANCE);
				if ( istar[obs] < 0 ) {
					errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
					mem_free(istar); return -1;
				}
				cell->LPcnt++;
			}
			else {
				/* (b) Use the argmax operator to identify the best basis */
				istar[obs] = computeIstar(prob->num, prob->coord, cell->basis, cell->sigma, cell->delta, Xvect, piCbarX,
						cell->omega->vals[obs]+offset, obs, cell->k, FALSE, &argmax, FALSE);
			}
		}
	}

	/* allocate memory to hold a new cut */
	cut = newCut(prob->num->prevCols, cell->omega->cnt, cell->omega->cnt);

	/* Go through all the cuts and form the coefficients using the basis indentified in the previous step. */
	for (obs = 0; obs < cell->omega->cnt; obs++) {
		for ( cnt = 0; cnt <= cell->basis->vals[istar[obs]]->phiLength; cnt++ ) {
			sigmaIdx  = cell->basis->vals[istar[obs]]->sigmaIdx[cnt];
			lambdaIdx = cell->basis->vals[istar[obs]]->lambdaIdx[cnt];
			if ( cnt == 0 )
					multiplier = 1.0;
				else
					multiplier = cell->omega->vals[obs][prob->num->rvbOmCnt+prob->num->rvCOmCnt+cell->basis->vals[istar[obs]]->omegaIdx[cnt]];

			cut->alpha += (cell->sigma->vals[sigmaIdx].pib + cell->delta->vals[lambdaIdx][obs].pib)*cell->omega->probs[obs]*multiplier;

			for (c = 1; c <= prob->num->cntCcols; c++)
				cut->beta[prob->coord->colsC[c]] += cell->sigma->vals[sigmaIdx].piC[c]*cell->omega->probs[obs]*multiplier;
			for (c = 1; c <= prob->num->rvCOmCnt; c++)
				cut->beta[prob->coord->rvCols[c]] += cell->delta->vals[lambdaIdx][obs].piC[c]*cell->omega->probs[obs]*multiplier;
		}
	}
	cut->alphaIncumb = cut->alpha;
	cut->beta[0] = 1.0;

	/* (c) add cut to the master problem  */
	if ( (cutIdx = addCut2Master(cell, cut, FALSE, prob->num->prevCols, 0.0)) < 0 ) {
		errMsg("algorithm", "formSDCut", "failed to add the new cut to master problem", 0);
		return -1;
	}

	mem_free(istar);
	return cutIdx;
}//END formCut()

BOOL optimalBenders(probType **prob, cellType *cell) {
	double gap;

	if ( cell->k > config.MIN_ITER ) {
		gap = (cell->ub - cell->candidEst)/abs(cell->candidEst);
		if ( gap <= config.EPSILON )
			return TRUE;
	}

#if defined(ALGO_CHECK)
	printf("Optimality gap = (%.6lf - %.6lf) = %.6lf\n", cell->ub, cell->candidEst, cell->ub - cell->candidEst);
#endif

	return FALSE;
}//END optimalBenders()

void writeBendersStatistic(FILE *soln, probType **prob, cellType *cell, string probName, int numStages) {
	int t;

	fprintf(soln, "\n\n====================================================================================================================================\n");
	fprintf(soln, "-------------------------------------------------------- Problem Information -------------------------------------------------------\n");
	fprintf(soln, "====================================================================================================================================\n");
	fprintf(soln, "Problem                            : %s\n", probName);
	fprintf(soln, "Number of stages                   : %d\n", numStages);
	for ( t = 0; t < numStages; t++ ) {
		fprintf(soln,  "------------------------------------------------------------------------------------------------------------------------------------\n");
		fprintf(soln,  "Stage %d\n", t);
		fprintf(soln,  "Number of decision variables (u_t) = %d\t\t", prob[t]->sp->mac);
		fprintf(soln,  "(Continuous = %d\tInteger = %d\tBinary = %d)\n", prob[t]->sp->mac - prob[t]->sp->numInt - prob[t]->sp->numBin, prob[t]->sp->numInt, prob[t]->sp->numBin);
		fprintf(soln,  "Number of constraints              = %d\n", prob[t]->sp->mar);
		if ( prob[t]->omegas != NULL ) {
			fprintf(soln,  "Number of random variables (omega) = %d\t\t", prob[t]->omegas->numRV);
			fprintf(soln,  "(a_t = %d; b_t = %d; c_t = %d; d_t = %d; A_t = %d; B_t = %d; C_t = %d; D_t = %d)\n", prob[t]->num->rvaOmCnt, prob[t]->num->rvbOmCnt, prob[t]->num->rvcOmCnt, prob[t]->num->rvdOmCnt,
					prob[t]->num->rvAOmCnt, prob[t]->num->rvBOmCnt, prob[t]->num->rvCOmCnt, prob[t]->num->rvDOmCnt);
		}
		else
			fprintf(soln,  "Number of random variables (omega) = 0\n");
	}

	fprintf(soln, "\n====================================================================================================================================\n");
	fprintf(soln, "----------------------------------------------------------- Optimization -----------------------------------------------------------\n");
	fprintf(soln, "====================================================================================================================================\n");
	fprintf(soln, "Algorithm                          : Benders Decomposition\n");
	fprintf(soln, "Number of iterations               : %d\n", cell->k);
	fprintf(soln, "Lower bound estimate               : %lf\n", cell->candidEst);
	fprintf(soln, "Upper bound estimate               : %lf\n", cell->ub);
	fprintf(soln, "Optimality gap estimate            : %lf (%.3lf%%)\n", cell->ub - cell->candidEst, 100*(cell->ub - cell->candidEst)/cell->candidEst);
	fclose(soln);

}//END WriteStat

