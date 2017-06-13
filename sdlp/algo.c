/*
  * algo.c
 *
 *  Created on: Apr 2, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "sdlp.h"

extern configType config;
extern string outputDir;

int algo(string probName, oneProblem *orig, stocType *stoc, timeType *tim) {
	probType **prob = NULL;
	cellType **cell = NULL;
	vector	 observ;
	FILE *solnFile;

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell) )
		goto TERMINATE;

	/* allocate memory to hold an observation */
	if ( !(observ = (vector) arr_alloc(stoc->numOmega, double)) )
		errMsg("allocation", "algo", "observ", 0);

	printAlgoDetails(0);
	while (TRUE) {
#if VERBOSE
		printf("\nIteration %4d :: ", cell[0]->k+1); fflush(stdout);
#else
		if ( (cell[0]->k) % 100 == 0)
			printf("\nIteration %4d :: ", cell[0]->k); fflush(stdout);
#endif

		/* if optimality conditions have been satisfied then break the while loop and exit. */
		if ( optimal(prob, cell, tim->numStages) )
			break;

		/* generate a new observation */
		generateOmega(stoc, observ, &config.RUN_SEED);

		/* forward pass */
		if (forwardPass(prob, cell, observ, tim->numStages)) {
			errMsg("algorithm", "algo","failed in forward pass", 0);
			goto TERMINATE;
		}

		/* backward pass */
		if (backwardPass(prob, cell, observ, tim->numStages)) {
			errMsg("algorithm", "algo", "failed in backward pass", 0);
			goto TERMINATE;
		}

		if ( config.QUADRATIC ) {
			/* check to see if there is an improvement */
			checkImprovement(prob, cell, tim->numStages);

		}
	}

	printf("\n\n-------------------------------------------------------------------------------------------------------------------------------\n");
	printf("Successfully completed excecution of SDLP algorithm on %s.\n", probName);
	printf("-------------------------------------------------------------------------------------------------------------------------------\n");

	/* Print solution details */
	solnFile = openFile(outputDir, "detailedSDLPsols.dat", "w");
	printSolutionDetails(solnFile, probName, prob, cell, tim->numStages);
	printSolutionShort(stdout, probName, prob, cell, tim->numStages);

	/* release memory allocated to different structures used in the algorithm */
	TERMINATE:
	cleanupAlgo(prob, cell, tim->numStages);
	mem_free(observ);

	return 0;
}//END algo()

int forwardPass(probType **prob, cellType **cell, vector observ, int numStages) {
	long long pathIdx = 0;
	int		t, status, obs, incumbIdx, pathOld = 0;

	/************************************************* setup and solve stage problems *****************************************************/
	/* since primal solution for terminal stage are not used we do not solve it here. Rather we solve it on the backward pass */
	for ( t = 0; t < numStages-1; t++ ) {
		cell[t]->k++;
#if VERBOSE
		printf("\nStage-%d :: ", t); fflush(stdout);
#endif

		/* update the right-hand side with state information for non-root stages */
		if ( t != 0 ) {
			/* update omega structure with the new observation */
			pathOld = cell[t]->omega->pathCurrent;
			obs = calcOmega(prob[t]->omegas, cell[t]->omega, observ+prob[t]->omegas->beg, pathIdx);
			pathIdx = cell[t]->omega->pathIdx[cell[t]->omega->pathCurrent];
			if ( pathOld != cell[t]->omega->pathCurrent )
				/* Current path is different from the path observed in previous iteration */
				cell[t]->incumb->chg = TRUE;

			/* copy the deterministic right-hand side and change it with endogenous state information */
			computeEndoRHS(prob[t]->bBar, prob[t]->Cbar, cell[t-1]->candidU, cell[t]->rhs);

			/* change the right-hand side with exogenous state information */
			if ( computeExoRHS(cell[t]->sp->lp, cell[t]->sda, prob[t]->coord, prob[t]->num, cell[t]->omega->vals[obs], cell[t-1]->candidU, cell[t]->rhs)) {
				errMsg("allocation", "forwardPass", "failed to change the right-hand side with uncertainty and state information", 0);
				return 1;
			}
		}
		else
			/* copy the deterministic right-hand side */
			computeEndoRHS(prob[t]->bBar, prob[t]->Cbar, NULL, cell[t]->rhs);

		/* select the incumbent solution to be used */
		incumbIdx = selectIncumb(cell[t]->incumb, cell[t]->omega);
#if VERBOSE
			printf("\tIncumbent chosen = %d.\t", incumbIdx); fflush(stdout);
#endif

		if ( config.QUADRATIC && cell[t]->incumb->chg ) {
			/* If decision simulation problem is solved as a quadratic program, then update the right-hand side and bounds using current incumbent
			 * solution */
			if ( changeQPrhs(cell[t]->sp->lp, prob[t+1]->coord->colsC, prob[t+1]->num->cntCcols, prob[t]->num->rows,
					prob[t]->Dbar, prob[t]->bBar, cell[t]->cuts, cell[t]->incumb->vals[incumbIdx], cell[t]->rhs, cell[t]->k, cell[t]->lb) ) {
				errMsg("algorithm", "forwardPass", "failed to change the proximal parameter", 0);
				return 1;
			}

			if ( changeQPbds(cell[t]->sp->lp, prob[t]->num->cols, prob[t]->sp->bdl, prob[t]->sp->bdu, cell[t]->incumb->vals[incumbIdx]) ) {
				errMsg("algorithm", "forwardPass", "failed to change the proximal parameter", 0);
				return 1;
			}

			/* update the proximal term for the subproblem */
			if ( constructQP(cell[t]->sp->lp, prob[t]->num->cols, cell[t]->incumb->quadScalar) ) {
				errMsg("algorithm", "forwardPass", "failed to change the proximal parameter", 0);
				return 1;
			}
#if VERBOSE
			printf("Update complete.\n");
#endif
		}
		else if (cell[t]->lbType == NONTRIVIAL ) {
			/* if the decision simulation problem is solved as a linear program and it has non-trivial lower bound, then update the right-hand side of
			 * minorants to reflect this. This update for when decision simulation problem is solved as quadratic program is performed in changeQPrhs */
			if ( updateCutsRHS(cell[t]->sp->lp, cell[t]->cuts, cell[t]->lb, cell[t+1]->k) ) {
				errMsg("algorithm", "forwardPass", "failed to change right-hand side to reflect non-trivial lower bound", 0);
				return 1;
			}
		}

		/* change coefficients of eta column */
		if ( changeEtaCol(cell[t]->sp->lp, prob[t]->num->cols, prob[t]->num->rows, cell[t]->k, cell[t]->cuts, cell[t]->lb) ) {
			errMsg("algorithm", "forwardPass", "failed to change the proximal parameter", 0);
			return 1;
		}

#ifdef ALGO_RUN
		char fname[NAMESIZE];
		sprintf(fname, "fProb%d_%d.lp", t, cell[t]->k);
		writeProblem(cell[t]->sp->lp, fname);
#endif

		/* solve the stage problem */
		if ( solveProblem(cell[t]->sp->lp, cell[t]->sp->name, cell[t]->sp->type, &status) ) {
			errMsg("solver", "forwardPass", "failed to solve stage problem", 0);
			return 1;
		}

		/* obtain the current primal solution */
		if ( getPrimal(cell[t]->sp->lp, cell[t]->candidU, prob[t]->num->cols) ) {
			errMsg("solver", "forwardPass", "failed to obtain primal solution for the stage problem", 0);
			return 1;
		}

		if ( config.QUADRATIC ) {
			/* find the norm of differences */
			cell[t]->incumb->normd_k = vXv(cell[t]->candidU, cell[t]->candidU, NULL, prob[t]->num->cols);
			if (cell[t]->k == 1)
				cell[t]->incumb->normd_k_1 = cell[t]->incumb->normd_k;

			/* Primal solution is \Delta u = u - \hat{u}, change it to u */
			addVectors(cell[t]->candidU, cell[t]->incumb->vals[incumbIdx], NULL, prob[t]->num->cols);

			/* Get the dual solution too */
			if ( getDual(cell[t]->sp->lp, cell[t]->pi, cell[t]->sp->mar+cell[t]->cuts->cnt) ) {
				errMsg("solver", "solveQPMaster", "failed to obtain dual solutions to master", 0);
				return 1;
			}
			if ( getDualSlacks(cell[t]->sp->lp, cell[t]->dj, prob[t]->num->cols) ) {
				errMsg("solver", "solveQPMaster", "failed to obtain dual slacks for master", 0);
				return 1;
			}
		}

		/* obtain the primal objective function value */
		cell[t]->candidEst = getObjective(cell[t]->sp->lp, PROB_LP);
	}

	return 0;
}//END forwardPass

int backwardPass(probType **prob, cellType **cell, vector observ, int numStages) {
	long long pathIdx;
	double	mubBar, futureVal;
	int 	t, extraRows = 0, idxSigma, idxCut;
	BOOL	newSigmaFlag;

	if (numStages == 2)
		pathIdx = 0;
	else
		pathIdx = cell[numStages-2]->omega->pathIdx[cell[numStages-2]->omega->pathCurrent];

	for ( t = numStages-1; t > 0; t-- ) {
		if ( t == numStages-1 ) {
#if VERBOSE
		printf("\nStage-%d :: ", t);
#endif
			/* update omega structure with the new observation, as this is not done in forward pass */
			cell[t]->k++;
			calcOmega(prob[t]->omegas, cell[t]->omega, observ+prob[t]->omegas->beg, pathIdx);

#if VERBOSE
			printf("\n-------------------------------------------------------------------------------------------------------------------------------\n");
#endif
			/* change the right-hand side with endogenous state */
			computeEndoRHS(prob[t]->bBar, prob[t]->Cbar, cell[t-1]->candidU, cell[t]->rhs);

			/* change the right-hand side with exogensous state */
			if ( computeExoRHS(cell[t]->sda, NULL, prob[t]->coord, prob[t]->num, cell[t]->omega->vals[cell[t]->omega->idx],
					cell[t-1]->candidU, cell[t]->rhs) ){
				errMsg("allocation", "backwardPass", "failed to change the right-hand side with uncertainty and state information", 0);
				return 1;
			}

			futureVal = 0.0;
		}
		else
			futureVal = cell[t]->cuts->vals[cell[t]->cuts->cnt-1]->alpha;

#ifdef ALGO_RUN
		char fname[NAMESIZE];
		sprintf(fname, "bProb%d_%d.lp", t, cell[t]->k);
		writeProblem(cell[t]->sda, fname);
#endif

		/* solve the stage dual approximation and obtain the dual solutions */
		if ( dualUpdates(cell[t]->sda, cell[t]->sp->name, prob[t]->num->rows+extraRows, prob[t]->num->cols, cell[t]->pi, &mubBar)) {
			errMsg("algorithm", "backwardPass","failed to complete d%ual updates", 0);
			return 1;
		}

		/* update all the stochastic components, indicate that the updates with respect to new node have been completed */
		idxSigma = stocUpdate(config.MAX_ITER, prob[t]->num, prob[t]->coord, prob[t]->Cbar, prob[t]->bBar, cell[t]->pi, mubBar, futureVal,
					cell[t]->lambda, cell[t]->sigma, &newSigmaFlag, cell[t]->delta, cell[t]->omega, cell[t]->k);

#ifdef STOC_CHECK
	double obj;
	obj = cell[t]->sigma->vals[idxSigma].pib - vXv(cell[t]->sigma->vals[idxSigma].piC, cell[t-1]->candidU, prob[t]->coord->colsC, prob[t]->num->cntCcols);
	obj += cell[t]->delta->vals[cell[t]->sigma->lambdaIdx[idxSigma]][cell[t]->omega->idx].pib - vXv(cell[t]->delta->vals[cell[t]->sigma->lambdaIdx[idxSigma]][cell[t]->omega->idx].piC,
			cell[t]->omega->vals[cell[t]->omega->idx], prob[t]->coord->rvCols, prob[t]->num->rvColCnt);
	printf("Objective function estimate at candidate solution = %lf\n", obj);
#endif

		/* form new optimality cut */
		idxCut = formCut(cell[t-1]->sp->lp, cell[t-1]->sda, cell[t], prob[t], cell[t-1]->cuts, cell[t-1]->candidU,
				prob[t-1]->num->rows, prob[t-1]->num->cols, t == (numStages - 1), numStages, cell[t-1]->pi, cell[t-1]->incumb->cutidx);
		if ( idxCut < 0 ) {
			errMsg("algorithm", "backwardPass", "failed to add the candidate cut", 0);
			return 1;
		}
		extraRows = 1;

		if ( numStages == 2 && (cell[t-1]->k % config.TAU == 0)) {
			formIncumbCut(cell[t], prob[t], cell[t-1]->sp->lp, cell[t-1]->sda, cell[t-1]->cuts, cell[t-1]->incumb->vals[cell[t-1]->incumb->idx],
					prob[t-1]->num->rows, prob[t-1]->num->cols, t == (numStages - 1), numStages, cell[t-1]->pi, cell[t-1]->incumb->cutidx);
		}
	}

	return 0;
}//END backwardPass()

void computeEndoRHS(sparseVector *bBar, sparseMatrix *Cbar, vector candidU, vector rhs) {
	int n;

	/* right-hand side: add the fixed part */
	for ( n = 1; n <= bBar->cnt; n++ )
		rhs[bBar->col[n]] = bBar->val[n];

	if ( candidU != NULL )
		/* transfer matrix: fixed part */
		rhs = MSparsexvSub(Cbar, candidU, rhs);

}//END computeEndoRHS()

int computeExoRHS(LPptr lp, LPptr sda, coordType *coord, numType *num, vector observ, vector candidut, vector rhs) {
	sparseVector bOmega;
	sparseMatrix COmega;
	intvec 		 indices;
	int 		 n, offset, status;

	if ( !(indices = (intvec) arr_alloc(num->rows, int)))
		errMsg("allocation", "computeRHS", "indices", 0);
	for ( n = 0; n < num->rows; n++) {
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

	/* change the right hand side in the solver for stage decision simulation problem */
	status = changeRHS(lp, num->rows, indices, rhs+1);
	if ( status ) {
		errMsg("solver", "computeRHS", "failed to change right-hand side in solver", 0);
		return 1;
	}

	if ( sda != NULL ) {
		/* change the right-hand side in the solver for the stage dual approximation problem */
		status = changeRHS(sda, num->rows, indices, rhs+1);
		if ( status ) {
			errMsg("solver", "computeRHS", "failed to change right-hand side in solver", 0);
			return 1;
		}
	}

	mem_free(indices);
	return 0;
}//END computeExoRHS()

int changeEtaCol(LPptr lp, int numCols, int numRows, int k, cutsType *cuts, double lb) {
	vector	coef;
	double	etaCoef[1], etaBds[1];
	int 	status, c, etaCol[1];
	char	bdsType[1];

	etaCol[0] = numCols;
	bdsType[0] = 'L';

	/* array of coefficients for the \eta column. */
	if (!(coef = (vector) arr_alloc(cuts->cnt, double)))
		errMsg("allocation", "chgEtaCol", "coef", 0);

	for (c = 0; c < cuts->cnt; c++) {
		/* Currently both incumbent and candidate cuts are treated similarly, and sunk as iterations proceed */
		coef[cuts->vals[c]->rowNum - numRows] = (double) (k-1) / (double) cuts->vals[c]->numObs;
	}

	/* change the eta column in the stage problem, which corresponds to the eta variable, starting with the row after the D matrix
	 * and ending with the row of the last cut. */
	if ( changeCol(lp, numCols, coef, numRows, numRows+ cuts->cnt) ) {
		errMsg("solver", "chgEtaCol", "failed to change eta column in the stage problem", 0);
		return 1;
	}

	/* if feasibility cut is added without any general cut, eta's lower bound should be changed to lower bound computed from mean value solution
	 * and its objective function coefficient should be zero. Once a general cut is encountered, revert to 1.0 for objective coefficient and -\infty
	 * for lower bound on eta */
	if ( cuts->cnt <= 1) {
		if ( cuts->cnt > 0 ) {
			etaCoef[0] = 1.0;
			etaBds[0]  = -INFBOUND;
		}
		else {
			etaCoef[0] = 0.0;
			etaBds[0] = lb;
		}
		status = changeObjx(lp, 1, etaCol, etaCoef);
		if ( status ) {
			errMsg("solver", "changeEtaCol", "failed to change the objective coefficient of eta column in objective function value", 0);
			return 1;
		}
		status = changeBDS(lp, 1, etaCol, bdsType, etaBds);
		if ( status ) {
			errMsg("solver", "changeEtaCol", "failed to change the bound for eta column", 0);
			return 1;
		}
	}

	mem_free(coef);
	return 0;
}//END chgEtaCol()

int updateCutsRHS(LPptr lp, cutsType *cuts, double lb, int numObs) {
	vector rhs;
	intvec indices;
	int cnt;

	if ( !(rhs = (vector) arr_alloc(cuts->cnt, double)) )
		errMsg("allocation", "updateRHS", "rhs", 0);
	if ( !(indices = (intvec) arr_alloc(cuts->cnt, int)) )
		errMsg("allocation", "updateRHS", "indices", 0);

	for ( cnt = 0; cnt < cuts->cnt; cnt++ ) {
		rhs[cnt] = cuts->vals[cnt]->alpha + ((double) numObs / (double) cuts->vals[cnt]->numObs - 1) * lb;
		indices[cnt] = cuts->vals[cnt]->rowNum;
	}

	if ( changeRHS(lp, cuts->cnt, indices, rhs) ) {
		errMsg("solver", "updateRHS", "failed to change the right-hand side", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(indices);

	return 0;
}//END updateCutsRHS

int dualUpdates(LPptr lp, string name, int numRows, int numCols, vector pi, double *mubBar) {
	int 	status;

	/* solve the terminal stage problem as a linear program */
	if ( solveProblem(lp, name, PROB_LP, &status) ) {
		errMsg("solver", "backPass", "failed to solve terminal stage problem", 0);
		return 1;
	}

#ifdef STOC_CHECK
	printf("Objective function value = %lf\t", getObjective(lp, PROB_LP)); fflush(stdout);
#endif

	/* obtain the dual solution */
	if (getDual(lp, pi, numRows) ) {
		errMsg("solver", "backwardPass", "failed to obtain optimal dual solutions to TDA problem", 0);
		return 1;
	}

	/* compute \bar{\mu} */
	if (computeMu(lp, numCols, mubBar) ) {
		errMsg("algorithm", "backwardPass", "failed to compute mu for stochastic updates", 0);
		return 1;
	}

	return 0;
}//END dualUpdates()

int computeMu(LPptr lp, int numCols, double *mubBar) {
	vector	dj, u;
	intvec	cstat;
	int		n;

	(*mubBar) = 0.0;

	if ( !(dj = (vector) arr_alloc(numCols+1, double)))
		errMsg("allocation", "computeMu", "dual slacks", 0);
	if ( !(u = (vector) arr_alloc(numCols+1, double)))
		errMsg("allocation", "computeMu", "TDA solutions", 0);

	if ( getPrimal(lp, u, numCols) ) {
		errMsg("solver", "forOptPass", "failed to obtain primal solution", 0);
		return 1;
	}
	if (getDualSlacks(lp, dj, numCols) ) {
		errMsg("solver", "computeMu", "failed to obtain dual slacks", 0);
		return 1;
	}

	/* extra column for eta if the stage problem is a QP */
	if ( !(cstat = (intvec) arr_alloc(numCols+2, int)) )
		errMsg("allocation", "computeMu", "column status", 0);
	if (getBasis(lp, cstat+1, NULL)) {
		errMsg("solver", "computeMu", "failed to get column status", 0);
		return 1;
	}

	for (n = 1; n <= numCols;  n++) {
		switch (cstat[n]) {
		case AT_LOWER:
			(*mubBar) += dj[n]*u[n];
			break;
		case AT_UPPER:
			(*mubBar) += dj[n]*u[n];
			break;
		default:
			break;
		}
	}

	mem_free(u); mem_free(cstat); mem_free(dj);

	return 0;
}//END computeMu()

void cleanupAlgo(probType **prob, cellType **cell, int T) {

	freeCellType(prob, cell, T);
	freeProbType(prob, T);

}//END cleanupAlgo()

void printAlgoDetails(int item) {

	printf("\n-------------------------------------------------------------------------------------------------------------------------------\n");
	printf("Starting stochastic dynamic linear programming algorithm\n");
	printf("Minimum iterations: %d\n", config.MIN_ITER);
	printf("-------------------------------------------------------------------------------------------------------------------------------\n");

}//END printAlgoDetails()

void printSolutionShort(void *fPtr, string probName, probType **prob, cellType **cell, int numStages) {

	fprintf(fPtr, "\n===============================================================================================================================\n");
	fprintf(fPtr, "Number of iterations                      = %d\n", cell[0]->k);
	fprintf(fPtr, "Objective function estimate at root stage = %lf\n", cell[0]->incumb->est[0]);
	fprintf(fPtr, "\n===============================================================================================================================\n");

}//END printSolutionShort()

void printSolutionDetails (void *fPtr, string probName, probType **prob, cellType **cell, int numStages) {
	int t, n;

	fprintf(fPtr, "\n=============================================================================================================\n");
	fprintf(fPtr, "Number of iterations                      = %d\n", cell[0]->k);
	fprintf(fPtr, "Objective function estimate at root stage = %lf\n", cell[0]->incumb->est[0]);

	for (t = 1; t < numStages; t++ ) {
		/* Details of stochastic elements */
		fprintf(fPtr, "-------------------------------------------------------------------------------------------------------------\n");
		fprintf(fPtr, "                                         --- Stage %d ---                                                    \n", t);
		fprintf(fPtr, "-------------------------------------------------------------------------------------------------------------\n");
		fprintf(fPtr, "Number of observations encountered = %d\n", cell[t]->omega->cnt);
		for ( n = 0; n < cell[t]->omega->cnt; n++)
			fprintf(fPtr, "%lf\t%lf\n", cell[t]->omega->vals[n][1] + prob[t]->omegas->mean[1], (double) cell[t]->omega->weights[n]/cell[t]->k);
		fprintf(fPtr, "\n");

		fprintf(fPtr, "Number of lambda's encountered = %d\n", cell[t]->lambda->cnt);
		for ( n = 0; n < cell[t]->lambda->cnt; n++){
			fprintf(fPtr, "%d: ", n);
			printVector(cell[t]->lambda->vals[n], prob[t]->num->rvRowCnt, fPtr);
		}
		fprintf(fPtr, "\n");

		fprintf(fPtr, "Number of sigma's encountered = %d\n", cell[t]->sigma->cnt);
		for ( n = 0; n < cell[t]->sigma->cnt; n++ ) {
			fprintf(fPtr, "%d: (%d,%d);\t%lf\t; ", n, cell[t]->sigma->ck[n], cell[t]->sigma->lambdaIdx[n], cell[t]->sigma->vals[n].pib);
			printVector(cell[t]->sigma->vals[n].piC, prob[t]->num->cntCcols, fPtr);
		}
		fprintf(fPtr, "\n");

		/* details of approximation */
		fprintf(fPtr, "Number of minorants in the approximation = %d\n", cell[t-1]->cuts->cnt);
		for ( n = 0; n < cell[t-1]->cuts->cnt; n++ ) {
			fprintf(fPtr, "(%d) ", cell[t-1]->cuts->vals[n]->numObs);
			printIntvec(cell[t-1]->cuts->vals[n]->iStar-1, cell[t-1]->cuts->vals[n]->numIstar, fPtr);
		}

		/* details of incumbent */
		if ( cell[t]->incumb != NULL ) {
			fprintf(fPtr, "Number of incumbents generated = %d\n", cell[t]->incumb->cnt-1);
			printVector(cell[t]->incumb->est, cell[t]->incumb->cnt-1, fPtr);
		}
	}
	fprintf(fPtr, "=============================================================================================================\n");

	fclose(fPtr);

}//END printSolutionDetails()
