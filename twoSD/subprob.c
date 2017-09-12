/*
 * subprob.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "stoc.h"

/* This function will solve a new subproblem. This involves replacing the right-hand side of the subproblem with new values, based upon some
 * observation of omega, and some X vector of primal variables from the master problem.  Generally, the latest observation is used.  When
 * forming a normal cut, the candidate x should be used, while the incumbent x should be used for updating the incumbent cut. */
int solveSubprob(probType *prob, oneProblem *subproblem, vector Xvect, basisType *basis, lambdaType *lambda, sigmaType *sigma, deltaType *delta, int deltaRowLength,
		omegaType *omega, int omegaIdx, BOOL newOmegaFlag, int currentIter, double TOLERANCE) {
	vector 	rhs, cost;
	intvec	indices;
	int  	status, n, offset = 0, basisIdx;

	if ( !(indices = (intvec) arr_alloc(max(prob->num->rows, prob->num->cols), int)) )
		errMsg("allocation", "solve_subporb", "indices", 0);
	for ( n = 0; n < max(prob->num->rows,prob->num->cols); n++ )
		indices[n] = n;

	/* (a) compute the right-hand side using current observation and first-stage solution */
	rhs = computeRHS(prob->num, prob->coord, prob->bBar, prob->Cbar, Xvect, omega->vals[omegaIdx]+offset);
	if ( rhs == NULL ) {
		errMsg("algorithm", "solveSubprob", "failed to compute subproblem right-hand side", 0);
		return 1;
	}

	/* (b) change the right-hand side in the solver */
	if ( changeRHS(subproblem->lp, prob->num->rows, indices, rhs + 1) ) {
		errMsg("solver", "solve_subprob", "failed to change the right-hand side in the solver",0);
		return 1;
	}

	offset = prob->num->rvbOmCnt + prob->num->rvCOmCnt;
	/* (c) compute the cost coefficients using current observation */
	cost = computeCostCoeff(prob->num, prob->coord, prob->dBar, omega->vals[omegaIdx], offset);
	if ( cost == NULL ) {
		errMsg("algorithm", "solveSubprob", "failed to compute subproblem cost coefficients", 0);
		return 1;
	}

	/* (d) change cost coefficients in the solver */
	if ( changeObjx(subproblem->lp, prob->num->cols, indices, cost+1) ) {
		errMsg("solver", "solve_subprob", "failed to change the cost coefficients in the solver",0);
		return 1;
	}

#if defined(ALGO_CHECK)
	writeProblem(subproblem->lp, "subproblem.lp");
#endif

	/* (e) Solve the subproblem to obtain the optimal dual solution. */
	if ( solveProblem(subproblem->lp, subproblem->name, subproblem->type, &status) ) {
		if ( status == STAT_INFEASIBLE ) {
			printf("Subproblem is infeasible: need to create feasibility cut.\n");
			return 1;
		}
		else {
			errMsg("algorithm", "solveSubprob", "failed to solve subproblem in solver", 0);
			return 1;
		}
	}

#if defined(STOCH_CHECK)
	double obj;
	obj = getObjective(subproblem->lp, PROB_LP);
	printf("Objective value of Subproblem  = %lf\n", obj);
#endif

	/* (f) update the stochastic elements in the problem */
	basisIdx = stochasticUpdates(prob, subproblem, basis, lambda, sigma, delta, deltaRowLength, omega, omegaIdx, newOmegaFlag, currentIter, TOLERANCE);
	if ( basisIdx < 0 ) {
		errMsg("algorithm", "solveSubprob", "stochastic updates failed", 0);
		return 1;
	}

#if defined(STOCH_CHECK)
	int sigmaIdx, lambdaIdx; double multiplier, obj1;
	obj1 = 0;
	for ( n = 0; n <= ->vals[basisIdx]->phiLength; n++ ) {
		sigmaIdx = ->vals[basisIdx]->sigmaIdx[n];
		lambdaIdx = ->vals[basisIdx]->lambdaIdx[n];
		if ( n == 0 )
			multiplier = 1.0;
		else
			multiplier = omega->vals[omegaIdx][prob->num->rvbOmCnt+prob->num->rvCOmCnt+->vals[basisIdx]->omegaIdx[n]];
		obj1 += multiplier*(cell->sigma->vals[sigmaIdx].pib - vXv(cell->sigma->vals[sigmaIdx].piC, Xvect, prob->coord->colsC, prob->num->cntCcols));
		obj1 += multiplier*(delta->vals[lambdaIdx][omegaIdx].pib - vXv(->vals[lambdaIdx][omegaIdx].piC,
				omega->vals[omegaIdx], prob->coord->rvCols, prob->num->rvCOmCnt));
	}
	printf("Objective function estimate    = %lf\n", obj1);
	if ( fabs(obj-obj1) > 0.001 )
		printf("WARNING: The objective function and the estimate computed using stochastic elements do not match.\n");
#endif

	mem_free(rhs); mem_free(cost); mem_free(indices);
	return 0;
}// END solveSubprob()

/* This function computes the right hand side of the subproblem, based on a given X vector and a given observation of omega.
 * It is defined as:
 * 			rhs = R(omega) - T(omega) x X
 * and is calculated as:
 * 			rhs = (Rbar - Tbar x X) + (Romega - Tomega x X)
 *
 * where the "bar" denotes the fixed or mean value, and the "omega" denotes a random variation from this mean. The function allocates an array
 * for the vector, which must be freed by the customer.  Also, the zeroth position of this rhs vector is reserved, and the actual values begin at rhs[1].
 * R is b, and T is C
 \***********************************************************************/
vector computeRHS(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, vector X, vector obs) {
	int cnt;
	vector rhs;
	sparseVector bomega;
	sparseMatrix Comega;

	bomega.cnt = num->rvbOmCnt;	bomega.col = coord->omegaRow; bomega.val=obs;

	Comega.cnt = num->rvCOmCnt; Comega.col = coord->omegaCol + num->rvbOmCnt;
	Comega.row = coord->omegaRow + num->rvbOmCnt; Comega.val = obs + num->rvbOmCnt;

	/* Start with the values of b(omega) -- both fixed and varying */
	rhs = expandVector(bBar->val, bBar->col, bBar->cnt, num->rows);
	for (cnt = 1; cnt <= bomega.cnt; cnt++)
		rhs[bomega.col[cnt]] += bomega.val[cnt];

	/* (cumulatively) subtract values of C(omega) x X -- both fixed and varying */
	rhs = MSparsexvSub(Cbar, X, rhs);
	rhs = MSparsexvSub(&Comega, X, rhs);

	return rhs;
}//END computeRHS()

vector computeCostCoeff(numType *num, coordType *coord, sparseVector *dBar, vector obs, int offset) {
	vector cost;
	sparseVector cOmega;
	int	cnt;

	cOmega.cnt = num->rvdOmCnt; cOmega.col = coord->omegaCol+offset; cOmega.val = obs+offset;

	cost = expandVector(dBar->val, dBar->col, dBar->cnt, num->cols);
	for (cnt = 1; cnt <= cOmega.cnt; cnt++)
		cost[cOmega.col[cnt]] += cOmega.val[cnt];

	return cost;
}//END computeCostCoeff()

void chgRHSwSoln(sparseVector *bBar, sparseMatrix *Cbar, vector rhs, vector X) {
	int cnt;

	/* copy the original right-hand side */
	for (cnt = 1; cnt <= bBar->cnt; cnt++)
		rhs[bBar->col[cnt]] = bBar->val[cnt];

	/* change the right-hand side with first stage solution */
	rhs = MSparsexvSub(Cbar, X, rhs);

}//END chgRHSwMean()

int chgRHSwObserv(LPptr lp, numType *num, coordType *coord, vector observ, vector spRHS, vector X) {
	sparseVector bomega;
	sparseMatrix Comega;
	vector 	rhs;
	intvec	indices;
	int		cnt, stat1;

	bomega.cnt = num->rvbOmCnt;	bomega.col = coord->omegaRow; bomega.val = observ;

	Comega.cnt = num->rvCOmCnt; Comega.col = coord->omegaCol + num->rvbOmCnt;
	Comega.row = coord->omegaRow + num->rvbOmCnt; Comega.val = observ + num->rvbOmCnt;

	if ( !(indices = (intvec) arr_alloc(num->rows, int)) )
		errMsg("allocation", "chgRHSwRand", "indices", 0);
	if ( !(rhs = (vector) arr_alloc(num->rows+1, double)) )
		errMsg("allocation", "chgRHSwRand", "rhs", 0);

	/* copy right-hand side modified with mean information */
	for ( cnt = 1; cnt <= num->rows; cnt++ ) {
		rhs[cnt] = spRHS[cnt];
		indices[cnt-1] = cnt-1;
	}

	/* change right-hand side with randomness in b */
	for (cnt = 1; cnt <= bomega.cnt; cnt++)
		rhs[bomega.col[cnt]] += bomega.val[cnt];

	/* change right-hand side with randomness in transfer matrix */
	rhs = MSparsexvSub(&Comega, X, rhs);

	/* change the right-hand side in the solver */
	stat1 = changeRHS(lp, num->rows, indices, rhs + 1);
	if ( stat1 ) {
		errMsg("solver", "chgRHSwRand", "failed to change the right-hand side in the solver",0);
		return 1;
	}

	mem_free(rhs); mem_free(indices);
	return 0;

}//END chgRHSwRand()

int chgObjxwObserv(LPptr lp, vector cost, intvec indices, int rvdOmCnt, vector observ) {
	vector vals;
	int n;

	if ( !(vals = (vector) arr_alloc(rvdOmCnt + 1, double)) )
		errMsg("allocation", "chgObjwObserv", "vals", 0);

	for ( n = 1; n <= rvdOmCnt; n++ )
		vals[n] = cost[n] + observ[n];

	if ( changeObjx(lp, rvdOmCnt, indices+1, vals+1) ) {
		errMsg("solver", "chgObjswObserv", "failed to change the cost coefficients in the solver",0);
		return 1;
	}

	mem_free(vals);
	return 0;
}//END chgObjwObserv()

oneProblem *newSubproblem(oneProblem *subprob) {

	/* since the basic structure of subproblem is not modified during the course of the algorithm, we just load it onto the solver */
	subprob->lp = setupProblem(subprob->name, subprob->type, subprob->mac, subprob->mar, subprob->objsen, subprob->objx, subprob->rhsx, subprob->senx,subprob->matbeg, subprob->matcnt, subprob->matind, subprob->matval, subprob->bdl, subprob->bdu, NULL, subprob->cname, subprob->rname, subprob->ctype);
	if ( subprob->lp == NULL ) {
		errMsg("Problem Setup", "new_subprob", "subprob",0);
		return NULL;
	}

#if defined(SETUP_CHECK)
	if (writeProblem(subprob->lp, "newSubproblem.lp") ){
		errMsg("solver", "newSubproblem", "failed to write subproblems to file", 0);
		return NULL;
	}
#endif

	return subprob;
}//END new_subprob
