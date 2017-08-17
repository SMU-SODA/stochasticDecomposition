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

#include "twoSD.h"

extern configType config;

/* This function will solve a new subproblem. This involves replacing the right-hand side of the subproblem with new values, based upon some
 * observation of omega, and some X vector of primal variables from the master problem.  Generally, the latest observation is used.  When
 * forming a normal cut, the candidate x should be used, while the incumbent x should be used for updating the incumbent cut. */
int solveSubprob(probType *prob, cellType *cell, vector Xvect, int omegaIdx, BOOL newOmegaFlag) {
	vector 	rhs, cost;
	intvec	indices;
	int  	status, n, offset = 0, basisIdx;
	BOOL	newBasisFlag;

    if ( !(indices = (intvec) arr_alloc(max(prob->num->rows, prob->num->cols), int)) )
        errMsg("allocation", "solve_subporb", "indices", 0);
    for ( n = 0; n < max(prob->num->rows,prob->num->cols); n++ )
    	indices[n] = n;

	/* (a) compute the right-hand side using current observation and first-stage solution */
    rhs = computeRHS(prob->num, prob->coord, prob->bBar, prob->Cbar, Xvect, cell->omega->vals[omegaIdx]+offset);
    if ( rhs == NULL ) {
        errMsg("algorithm", "solveSubprob", "failed to compute subproblem right-hand side", 0);
        return 1;
    }

    /* (b) change the right-hand side in the solver */
    if ( changeRHS(cell->subprob->lp, prob->num->rows, indices, rhs + 1) ) {
        errMsg("solver", "solve_subprob", "failed to change the right-hand side in the solver",0);
        return 1;
    }

    offset = prob->num->rvbOmCnt + prob->num->rvCOmCnt;
    /* (c) compute the cost coefficients using current observation */
    cost = computeCostCoeff(prob->num, prob->coord, prob->cBar, cell->omega->vals[omegaIdx], offset);
    if ( cost == NULL ) {
    	errMsg("algorithm", "solveSubprob", "failed to compute subproblem cost coefficients", 0);
    	return 1;
    }

    /* (d) change cost coefficients in the solver */
    if ( changeObjx(cell->subprob->lp, prob->num->cols, indices, cost+1) ) {
    	errMsg("solver", "solve_subprob", "failed to change the cost coefficients in the solver",0);
    	return 1;
    }

#if defined(ALGO_CHECK)
    writeProblem(cell->subprob->lp, "subproblem.lp");
#endif

    /* (e) Solve the subproblem to obtain the optimal dual solution. */
    if ( solveProblem(cell->subprob->lp, cell->subprob->name, cell->subprob->type, &status) ) {
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
    obj = getObjective(cell->subprob->lp, PROB_LP);
    printf("Objective value of Subproblem  = %lf\n", obj);
#endif

    if ( getDual(cell->subprob->lp, cell->piS, prob->num->rows) ) {
        errMsg("algorithm", "solveSubprob", "failed to get the dual", 0);
        return 1;
    }
    if ( computeMU(cell->subprob->lp, prob->num->cols, &cell->mubBar) ) {
        errMsg("algorithm", "solveSubprob", "failed to compute mubBar for subproblem", 0);
        return 1;
    }
    basisIdx = calcBasis(cell->subprob->lp, prob->num->cols, prob->num->rows, cell->basis, &newBasisFlag);

	/* (f) update the stochastic elements in the problem */
	status = stochasticUpdates(prob->num, prob->coord, prob->bBar, prob->Cbar, cell->lambda, cell->sigma,
			cell->delta, cell->omega, newOmegaFlag, omegaIdx, config.MAX_ITER, cell->k, cell->piS, cell->mubBar);

#if defined(STOCH_CHECK)
	obj = cell->sigma->vals[status].pib - vXv(cell->sigma->vals[status].piC, Xvect, prob->coord->colsC, prob->num->cntCcols);
	obj += cell->delta->vals[cell->sigma->lambdaIdx[status]][omegaIdx].pib - vXv(cell->delta->vals[cell->sigma->lambdaIdx[status]][omegaIdx].piC,
			cell->omega->vals[omegaIdx], prob->coord->rvCols, prob->num->rvColCnt);
	printf("Objective function estimate    = %lf\n", obj);
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

    if (!(rhs =(vector) arr_alloc(num->rows+1, double)))
        errMsg("allocation", "computeRhs", "rhs",0);

    /* Start with the values of b(omega) -- both fixed and varying */
    for (cnt = 1; cnt <= bBar->cnt; cnt++)
        rhs[bBar->col[cnt]] +=  bBar->val[cnt];
    for (cnt = 1; cnt <= bomega.cnt; cnt++)
        rhs[bomega.col[cnt]] += bomega.val[cnt];

    /* (cumulatively) subtract values of C(omega) x X -- both fixed and varying */
    rhs = MSparsexvSub(Cbar, X, rhs);
    rhs = MSparsexvSub(&Comega, X, rhs);

    return rhs;
}//END computeRHS()

vector computeCostCoeff(numType *num, coordType *coord, sparseVector *cBar, vector obs, int offset) {
	vector cost;
	sparseVector cOmega;
	int	cnt;

	cOmega.cnt = num->rvcOmCnt; cOmega.col = coord->omegaCol+offset; cOmega.val = obs+offset;

	if ( !(cost = (vector) arr_alloc(num->cols+1, double)) )
		errMsg("allocation", "computeCostCoeff", "cost", 0);

	for (cnt = 1; cnt <= cBar->cnt; cnt++)
		cost[cBar->col[cnt]] = cBar->val[cnt];
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
