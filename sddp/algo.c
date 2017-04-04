/*
 * algo.c
 *
 *  Created on: Mar 27, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *     Contact: harsha@smu.edu
 *
 */

#include "sddp.h"
extern configType config;

int algo (oneProblem *orig, stocType *stoc, timeType *tim) {
	probType **prob = NULL;
	cellType **cell = NULL;

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell) )
		goto TERMINATE;

	printAlgoDetails(0);
	while (TRUE) {
		/* if optimality conditions have been satisfied then break the while loop and exit. */
		if ( optimal(prob, cell, tim->numStages) )
			break;

		/* forward pass */
		if (forwardPass(prob, cell, tim->numStages)) {
			errMsg("algorithm", "algo","failed in forward pass", 0);
			goto TERMINATE;
		}

		/* backward pass */
		if (backwardPass(prob, cell, tim->numStages)) {
			errMsg("algorithm", "algo", "failed in backward pass", 0);
			goto TERMINATE;
		}
	}

	/* release memory allocated to different structures used in the algorithm */
	TERMINATE:
	cleanupAlgo(prob, cell, tim->numStages);

	return 0;
}//END algo()

int forwardPass(probType **prob, cellType **cell, int numStages) {
	int		t, status, stat1, obs;

	/************************************************* setup and solve stage problems *****************************************************/
	/* since primal solution for terminal stage are not used we do not solve it here. Rather we solve it on the backward pass */
	for ( t = 0; t < numStages-1; t++ ) {
		cell[t]->k++;

		/* update the right-hand side with state information for non-root stages */
		if ( t != 0 ) {
			/* generate omega */
			obs = randInteger(&config.FORWPASS_SEED, cell[t]->omega->cnt);

			/* change the right-hand side with endogenous state information */
			computeEndoRHS(prob[t]->bBar, prob[t]->Cbar, cell[t-1]->candidU, cell[t]->rhs);

			/* change the right-hand side with exogenous state information */
			status = computeExoRHS(cell[t]->sp->lp, prob[t]->coord, prob[t]->num, cell[t]->omega->vals[obs], cell[t-1]->candidU, cell[t]->rhs);
			if ( status ) {
				errMsg("allocation", "forwardPass", "failed to change the right-hand side with uncertainty and state information", 0);
				return 1;
			}
		}

#ifdef ALGO_RUN
		char fname[NAMESIZE];
		sprintf(fname, "fProb%d_%d.lp", t, cell[t]->k);
		writeProblem(cell[t]->sp->lp, fname);
#endif

		/* solve the stage problem */
		status = solveProblem(cell[t]->sp->lp, cell[t]->sp->name, PROB_LP, &stat1);
		if (status) {
			errMsg("solver", "forwardPass", "failed to solve stage problem", 0);
			return 1;
		}

		/* obtain the current primal solution */
		stat1 = getPrimal(cell[t]->sp->lp, cell[t]->candidU, prob[t]->num->cols);
		if ( stat1 ) {
			errMsg("solver", "forwardPass", "failed to obtain primal solution for the stage problem", 0);
			return 1;
		}

		/* obtain the primal objective function value */
		cell[t]->candidEst = getObjective(cell[t]->sp->lp, PROB_LP);
	}

	return 0;
}//END forwardPass

int backwardPass(probType **prob, cellType **cell, int numStages) {
	int 	t, obs, status, stat1;
	intvec	iStar;

	/* solve the terminal problem as it was not solved during forward pass */
	t = numStages - 1;
	cell[t]->k++;
	computeEndoRHS(prob[t]->bBar, prob[t]->Cbar, cell[t-1]->candidU, cell[t]->rhs);

	for ( t = numStages-1; t > 0; t-- ) {
		if ( !(iStar = (intvec) arr_alloc(cell[t]->omega->cnt, int)) )
			errMsg("allocation", "backwardPass", "iStar", 0);

		for ( obs = 0; obs < cell[t]->omega->cnt; obs++ ) {
			/* select the subproblems to solve for cut generation */
			if ( randUniform(&config.BACKPASS_SEED) <= config.BACKPASS_PCT ) {
				/* change the right-hand side with exogenous information */
				status = computeExoRHS(cell[t]->sp->lp, prob[t]->coord, prob[t]->num, cell[t]->omega->vals[obs], cell[t-1]->candidU, cell[t]->rhs);
				if ( status ) {
					errMsg("allocation", "backwardPass", "failed to change the right-hand side with uncertainty and state information", 0);
					return 1;
				}

#ifdef ALGO_RUN
				char fname[NAMESIZE];
				sprintf(fname, "rProb%d_%d.lp", t, cell[t]->k);
				writeProblem(cell[t]->sp->lp, fname);
#endif

				/* solve the problem */
				status = solveProblem(cell[t]->sp->lp, cell[t]->sp->name, PROB_LP, &stat1);
				if (status) {
					errMsg("solver", "backwardPass", "failed to solve stage problem", 0);
					return 1;
				}

				/* update the stochastic elements */
				iStar[obs] = stocUpdate(prob[t], cell[t], cell[t-1]->candidU, obs);
				if ( iStar[obs] < 0 ) {
					errMsg("algorithm", "backwardPass", "failed to update stochastic elements", 0);
					return 1;
				}
			}
			else {
				/* mark the observation for which the best lower bound has to be  found */
				iStar[obs] = -1;
			}
		}

		/* identify the best lower bound for observations which were not solved */
		for ( obs = 0; obs < cell[t]->omega->cnt; obs++ )
			if ( iStar[obs] == -1 )
				iStar[obs] = computeIstar(prob[t]->num,prob[t]->coord, cell[t]->lambda, cell[t]->sigma, cell[t]->omega, cell[t]->delta,
						obs, cell[t-1]->candidU);

		/* create new optimality cut */
		stat1 = cell[t-1]->cuts->cnt;
		status = formOptCut(prob[t], cell[t], iStar, cell[t-1]->sp->lp, prob[t-1]->num->rows, prob[t-1]->num->cols, cell[t-1]->cuts, cell[t-1]->candidU);
		if ( status > stat1 )
			cell[t-1]->sp->mar++;

#ifdef CUT_CHECK
		double obj;
		obj = cell[t-1]->cuts->vals[cell[t-1]->cuts->cnt-1]->alpha + vXv(cell[t-1]->cuts->vals[cell[t-1]->cuts->cnt-1]->beta,
				cell[t-1]->candidU, prob[t]->coord->colsC, prob[t]->num->cntCcols);
		printf("Cut estimate = %lf\n", obj);
#endif

		/* clean up the stochastic elements */
		freeLambdaType(cell[t]->lambda, FALSE);
		freeSigmaType(cell[t]->sigma, FALSE);
		freeDeltaType(cell[t]->delta, cell[t]->omega->cnt, FALSE);
	}

	return 0;
}//END backwardPass

void computeEndoRHS(sparseVector *bBar, sparseMatrix *Cbar, vector candidU, vector rhs) {
	int n;

	/* right-hand side: add the fixed part */
	for ( n = 1; n <= bBar->cnt; n++ )
		rhs[bBar->col[n]] = bBar->val[n];

	/* transfer matrix: fixed part */
	rhs = MSparsexvSub(Cbar, candidU, rhs);

}//END computeEndoRHS()

int computeExoRHS(LPptr lp, coordType *coord, numType *num, vector observ, vector candidut, vector endoRHS) {
	sparseVector bOmega;
	sparseMatrix COmega;
	vector		 rhs;
	intvec 		 indices;
	int 		 n, offset, status;

	if ( !(indices = (intvec) arr_alloc(num->rows, int)))
		errMsg("allocation", "computeRHS", "indices", 0);
	if ( !(rhs = (vector) arr_alloc(num->rows+1, double)))
		errMsg("allocation", "computeRHS", "rhs", 0);
	for ( n = 0; n < num->rows; n++) {
		rhs[n+1] = endoRHS[n+1];
		indices[n]= n;
	}

	/* initialize bOmega and COmega; It is assumed that the realization vector is ordered as right-hand side, transfer matrix, cost coefficients */
	bOmega.cnt = num->rvbOmCnt;		bOmega.col = coord->omegaRow;			bOmega.val = observ;
	offset = num->rvbOmCnt;
	COmega.cnt = num->rvCOmCnt;		COmega.col = coord->omegaCol+offset;	COmega.row = coord->omegaRow+offset;	COmega.val = observ+offset;

	/* right-hand side: exogenous information */
	addVectors(rhs, bOmega.val, bOmega.col, bOmega.cnt);

	/* randomness in transfer matrix */
	rhs = MSparsexvSub(&COmega, candidut, rhs);

	/* change the right hand side in the solver for stage problem */
	status = changeRHS(lp, num->rows, indices, rhs+1);
	if ( status ) {
		errMsg("solver", "computeRHS", "failed to change right-hand side in solver", 0);
		return 1;
	}

	mem_free(rhs); mem_free(indices);
	return 0;
}//END computeExoRHS()

void printAlgoDetails(int item) {

	printf("\n-------------------------------------------------------------------------------------------------------------\n");
	printf("Starting stochastic dual dynamic programming algorithm\n");
	printf("Minimum iterations: %d\n", config.MIN_ITER);
	printf("-------------------------------------------------------------------------------------------------------------\n");

}//END printAlgoDetails()
