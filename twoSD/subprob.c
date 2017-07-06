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


//
//  subprob.c
//  multiAgentSP
//
//  Created by Shasha Wang on 2/4/16.
//  Copyright © 2016 Shasha Wang. All rights reserved.
//

#include "utils.h"
#include "smps.h"
#include "solver.h"
#include "prob.h"
#include "multiAgentSP.h"

extern configType config;

/* This function will solve a new subproblem. This involves replacing the right-hand side of the subproblem with new values, based upon some
 * observation of omega, and some X vector of primal variables from the master problem.  Generally, the latest observation is used.  When
 * forming a normal cut, the candidate x should be used, while the incumbent x should be used for updating the incumbent cut. */
int solveSubprob(probType *prob, cellType *cell, vector Xvect, vector observ) {
    intvec 	ind;
    vector 	rhs;
    int 	k, stat1, stat2;
    int status;
    clock_t tic, toc;

    if ( !(ind = (intvec) arr_alloc(prob->num->rows, int)) )
        errMsg("allocation", "solve_subporb", "ind", 0);

    for(k = 0; k < prob->num->rows; k++)
        ind[k] = k;

    /* compute the right-hand side using current observation and first-stage solution */
    rhs = computeRHS(prob->num, prob->coord, prob->bBar, prob->Cbar, Xvect, observ, cell->ID);
    if ( rhs == NULL ) {
        errMsg("algorithm", "solveSubprob", "failed to compute subproblem right-hand side", 0);
        return 1;
    }

#ifdef RHS_CHECK
    printf("\t\t\t    Subproblem RHS: ");
    printVector(rhs, prob->num->rows, NULL);
    writeProblem(cell->sp->lp, "subprobBefChgRHS.lp");
#endif

    /* change the right-hand side in the solver */
    stat1 = changeRHS(cell->sp->lp, prob->num->rows, ind, rhs + 1);
    if ( stat1 ) {
        errMsg("solver", "solve_subprob", "failed to change the right-hand side in the solver",0);
        return 1;
    }

#ifdef RHS_CHECK
    writeProblem(cell->sp->lp, "subprobAftChgRHS.lp");
#endif

    tic = clock();

    stat1 = solveProblem(cell->sp->lp, cell->sp->name, cell->sp->type, &stat2);

    toc = clock();
    cell->runTime->iterSolTime = ((double) (toc-tic)) / CLOCKS_PER_SEC;

    if ( stat1 ) {
        if ( stat2 == STAT_INFEASIBLE ) {
            printf("Subproblem is infeasible: need to create feasibility cut.\n");
            printf("cell = %d\n", cell->ID);
            cell->feasFlag = FALSE;
            return 0;
        }
        else {
            errMsg("algorithm", "solveSubprob", "failed to solve subproblem in solver", 0);
            return 1;
        }
    }
    else
        cell->feasFlag = TRUE;
    cell->LPcnt++;

#ifdef STOCH_CHECK
    double objV;
    objV = getObjective(cell->sp->lp, PROB_LP);
    printf("\t\t\t    Objective value of Subproblem-%d  = %lf\n", cell->ID, objV);
#endif

#ifdef SUB_SOLUTION
    /* record the objective function value */
    cell->optValM = getObjective(cell->sp->lp, PROB_LP);
    printf("\n%lf\n",cell->optValM);
#endif

    status=getDual(cell->sp->lp, cell->pi, prob->num->rows);
    if(status) {
        errMsg("algorithm", "solveSubprob", "failed to get the dual", 0);
        return 1;
    }

#ifdef STOCH_CHECK
    printf("\t\t\t    Dual: ");
    printVector(cell->pi, prob->num->rows, NULL);
#endif
    /* mubBar used in stochastic updates */
    stat1 = computeMU(cell->sp->lp, prob->num->cols, &cell->mubBar);
    if ( stat1 ) {
        errMsg("algorithm", "solveSubprob", "failed to compute mubBar for subproblem", 0);
        return 1;
    }

    mem_free(rhs);
    mem_free(ind);

    return 0;
}//END solve_subprob

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
vector computeRHS(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, vector X, vector obs, int ID) {
    int cnt;
    vector rhs;
    sparseVector bomega;
    sparseMatrix Comega;

    bomega.cnt = num->rvbOmCnt;	bomega.col = coord->omegaRow; bomega.val=obs;

    Comega.cnt = num->rvCOmCnt; Comega.col = coord->omegaCol + num->rvbOmCnt;
    Comega.row = coord->omegaRow + num->rvbOmCnt; Comega.val = obs + num->rvbOmCnt;

    if (!(rhs =(vector) arr_alloc(num->rows+1, double)))
        errMsg("Allocation", "computeRhs", "rhs",0);

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

void chgRHSwMean(sparseVector *bBar, sparseMatrix *Cbar, vector rhs, vector X) {
    int cnt;

    /* copy the original right-hand side */
    for (cnt = 1; cnt <= bBar->cnt; cnt++)
        rhs[bBar->col[cnt]] = bBar->val[cnt];

    /* change the right-hand side with first stage solution */
    rhs = MSparsexvSub(Cbar, X, rhs);

}//END chgRHSwMean()

int chgRHSwRand(LPptr lp, numType *num, coordType *coord, vector observ, vector spRHS, vector X, int agent) {
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

#if 0
    if (agent ==2)
        writeProblem(lp, "chgrhsevl.lp");
#endif

    /* change the right-hand side in the solver */
    stat1 = changeRHS(lp, num->rows, indices, rhs + 1);
    if ( stat1 ) {
        errMsg("solver", "chgRHSwRand", "failed to change the right-hand side in the solver",0);
        return 1;
    }

#if 0
    if (agent ==2)
        writeProblem(lp, "chgrhsevl1.lp");
#endif

    mem_free(rhs); mem_free(indices);
    return 0;

}//END chgRHSwRand()

cellType *newSubprob(probType *subprob, int agent, double weight) {
    cellType    *scell;
    int         length, i;

    if (!(scell = (cellType *) mem_malloc (sizeof(cellType))))
        errMsg("Memory allocation", "new_master", "Faile to allocate memory to master cell", 0);
    scell->ID = agent;
    scell->weight = weight;

    /* since the basic structure of subproblem is not modified during the course of the algorithm, we just load it onto the solver */
    subprob->sp->lp = setupProblem(subprob->sp->name, subprob->sp->type, subprob->sp->mac, subprob->sp->mar, subprob->sp->objsen, subprob->sp->objx, subprob->sp->rhsx, subprob->sp->senx,subprob->sp->matbeg, subprob->sp->matcnt, subprob->sp->matind, subprob->sp->matval, subprob->sp->bdl, subprob->sp->bdu, NULL, subprob->sp->cname, subprob->sp->rname, subprob->sp->ctype);
    if ( subprob->sp->lp == NULL ) {
        errMsg("Problem Setup", "new_subprob", "subprob->sp",0);
        return NULL;
    }

    scell->sp = subprob->sp;

#if 0
    int     status;
    char probName[NAMESIZE];
    sprintf(probName,"newSubprob%d.lp", agent);
    status = writeProblem(scell->sp->lp, probName);
    if ( status ) {
        errMsg("write problem", "new_subprob", "failed to write subproblems problem to file",0);
        return NULL;
    }
#endif

    scell->maxCuts = config.CUT_MULT * MAX_CUTS(subprob->num->prevCols);

    length = config.MAX_ITERATION + config.MAX_ITERATION / config.TAU + 1;


    if (!(scell->candidX = (vector) arr_alloc (subprob->num->cols+1,double)))
        errMsg("allocation", "duplicArray", "b", 1);

    scell->candidEst = 0.0;

    scell->dualStableFlag = FALSE;

    /* incumbent cuts are dropped only after adding to this structure */
    scell->cuts = newCuts(scell->maxCuts + 2);
    scell->fcuts = NULL;
    scell->incumbEst = 0.0;
    scell->quadScalar = 0.0;
    scell->LPcnt = 0;
    scell->feasCnt = 0;
    scell->infeasIncumb = FALSE;
    scell->feasFlag = TRUE;
    scell->incumbX = NULL;
    scell->incumbStdev = 0.0;
    scell->incumbChg = FALSE;
    scell->iCutIdx = 0;
    scell->iCutUpdt = 0;
    scell->gamma = 0.0;
    scell->normDk_1 = 0.0;
    scell->normDk = 0.0;
    scell->optFlag = FALSE;
    scell->k = 0;
    scell->newOmegaFlag = FALSE;
    scell->optValM =0.0;
    //MARK: check
    scell->pushcnt = 0;

    if ( !(scell->spRHS = (vector) arr_alloc(subprob->bBar->cnt + 1, double)) )
        errMsg("allocation", "newSoln", "scell->pi", 0);

    if ( !(scell->pi = (vector) arr_alloc(subprob->num->rows + 1, double)) )
        errMsg("allocation", "newSoln", "scell->pi", 0);

    if ( !(scell->pi_ratio = (vector) arr_alloc(config.SCAN_LEN, double)) )
        errMsg("allocation", "newMaster", "mcell->pi_ratio", 0);

    for (i = 0; i < config.SCAN_LEN; i++)
        scell->pi_ratio[i] = 0;

    if ( !(scell->di = (vector) arr_alloc(subprob->num->cols + 2, double)) ) //1 for eta, 1 for summation
        errMsg("allocation", "newSoln", "scell->di", 0);

    scell->lambda = newLambda(length, 0, subprob->num->rvRowCnt);
    scell->sigma = newSigma(length, subprob->num->rvColCnt, 0);
    scell->delta = newDelta(length);
    scell->omega = newOmega(config.MAX_ITERATION);

    /* initialize the time structure */
    if (!(scell->runTime =(runTimeType *) mem_malloc(sizeof(runTimeType))))
        errMsg("Allocation", "newMaster", "runTime", 0);
    scell->runTime->iterCutGenTime = 0.0; scell->runTime->iterSolTime = 0.0; scell->runTime->iterTime = 0.0;
    scell->runTime->totCutGenTime = 0.0; scell->runTime->totSolTime = 0.0; scell->runTime->totTime = 0.0;


    return scell;
}//END new_subprob
