/*
 * master.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */


///
//  master.c
//  multiAgentSP
//
//  Created by Shasha Wang on 1/27/16.
//  Copyright © 2016 Shasha Wang. All rights reserved.
//

#include "utils.h"
#include "smps.h"
#include "solver.h"
#include "prob.h"
#include "multiAgentSP.h"

extern configType config;
extern int numAgents;

/* This function is the regularized QP version of master problem. The master problem is solved after the newest cut is added to master problem,
 the incumbent cut is updated if necessary. Here the coefficients on all the cuts are updated, and finally master problem is solved. */
int solveQPMaster(numType *num, sparseVector *dBar, cellType **cell, int IniRow, double lb) {
	double 	d2 = 0.0; /* height at the candidate solution. */
	int 	status, stat1, i;
	clock_t tic, toc;

	if ( config.MULTI_CUT ) {
		for (i = 1; i < numAgents; i++) {
			status = changeEtaCol(cell[0]->sp->lp, num->rows, num->cols, cell[0]->k, cell[i]->cuts, lb, i, cell[i]->weight);
			if ( status ) {
				errMsg("algorithm", "solveMaster", "failed to change the eta column coefficients", 0);
				return 1;
			}
		}
	}
	else {
        //MARK: not sure we should pass 0 or 1
		status = changeEtaCol(cell[0]->sp->lp, num->rows, num->cols, cell[0]->k, cell[0]->cuts, lb, 1, cell[0]->weight);
		if ( status ) {
			errMsg("algorithm", "solveMaster", "failed to change the eta column coefficients", 0);
			return 1;
		}
	}

	if ( cell[0]->lbType == NONTRIVIAL ) {
		/* update the right-hand side of cuts to reflect the non-trivial lower bound */
		status = updateRHS(cell[0], IniRow, cell[0]->k);
		if ( status ) {
			errMsg("algorithm", "solveQPMaster", "failed to update right-hand side with lower bound information", 0);
			return 1;
		}
	}

	if ( cell[0]->incumbChg )
		cell[0]->incumbChg = FALSE;

#if 0
	writeProblem(cell[0]->sp->lp, "masterCell.lp");
#endif

	/* solve the master problem */
	tic = clock();

	status = solveProblem(cell[0]->sp->lp, cell[0]->sp->name, config.MASTERTYPE, &stat1);

	toc = clock();
    cell[0]->runTime->iterSolTime = ((double) (toc-tic)) / CLOCKS_PER_SEC;

	if ( status ) {
		writeProblem(cell[0]->sp->lp, "error.lp");
		errMsg("algorithm", "solveMaster", "failed to solve the master problem", 0);
		return 1;
	}

	/* increment the number of problems solved during algorithm */
	cell[0]->LPcnt++;

	/* record the objective function value */
	cell[0]->optValM = getObjective(cell[0]->sp->lp, cell[0]->sp->type);

	/* Get the most recent optimal solution to master program */
	status = getPrimal(cell[0]->sp->lp, cell[0]->candidX, num->cols);
	if ( status ) {
		errMsg("algorithm", "solveMaster", "failed to obtain the primal solution for master", 0);
		return 1;
	}

	/* add the incumbent back to change from \Delta X to X */
	for (i = 1; i <= num->cols; i++)
		d2 += cell[0]->candidX[i] * cell[0]->candidX[i];
	addVectors(cell[0]->candidX, cell[0]->incumbX, NULL, num->cols);

	/* update d_norm_k in soln_type. */
	if (cell[0]->k == 1)
		cell[0]->normDk_1 = d2;
	cell[0]->normDk = d2;

	/* Get the dual solution too */
	status = getDual(cell[0]->sp->lp, cell[0]->pi, cell[0]->sp->mar);
	if ( status ) {
		errMsg("solver", "solveQPMaster", "failed to obtain dual solutions to master", 0);
		return 1;
	}
	status = getDualSlacks(cell[0]->sp->lp, cell[0]->di, num->cols);
	if ( status ) {
		errMsg("solver", "solveQPMaster", "failed to obtain dual slacks for master", 0);
		return 1;
	}

	/* Find the highest cut at the candidate solution. where cut_height = alpha - beta(xbar + \Delta X) */
	cell[0]->candidEst = vXvSparse(cell[0]->candidX, dBar);

	if (cell[0]->sp->mar - IniRow > 0) {
		for ( i = 1; i < numAgents; i++){
			cell[0]->candidEst += cell[i]->weight * maxCutHeight(cell[0]->lbType, cell[i]->cuts, cell[0]->k, cell[0]->candidX, num->cols, lb, &status);
		}
	}
	else {
		if (cell[0]->lbType == TRIVIAL)
			cell[0]->candidEst += 0.0;
		else
			cell[0]->candidEst += lb;
	}


	/* Calculate gamma for next improvement check on incumbent x. */
	cell[0]->gamma = cell[0]->candidEst - cell[0]->incumbEst;

	return 0;
}//END solveQPMaster()

int addCuts2Master(cellType **cell, int lenX, double lb) {
    oneCut	*cut;
	intvec 	indices;
	int 	cnt, i, status, c;
    BOOL    dropIncumb = FALSE;

	if (!(indices = arr_alloc(lenX + 1, int)))
		errMsg("Allocation", "addcut2Master", "fail to allocate memory to coefficients of beta",0);
	for (cnt = 1; cnt <= lenX; cnt++)
		indices[cnt] = cnt - 1;

	if (config.MULTI_CUT) {
		/* one cut is added for every agent is added to the master */
        for ( i = 1; i < numAgents; i++ ) {
			/* identify the column index for eta corresponding to agent-i */
			indices[0] = lenX + cell[i]->ID - 1;

			/* if an incumbent cut was created in the current iteration, then two cuts need to be added to the master. */
			if ((cell[0]->k - cell[i]->iCutUpdt) % config.TAU == 0) {
                //cell[i]->pushcnt++;

				/* The incumbent cut is added first. Check to see if there is room for the incumbent cut, else drop the previous
				 * incumbent cut. The incumbent cut id does not change in this case. */
				if (cell[i]->cuts->cnt > cell[i]->maxCuts) {
                //if (cell[i]->pushcnt > cell[i]->maxCuts) {
					status = dropCut(cell, cell[i]->iCutIdx, i, cell[0]->sp->lp);
					if ( status ) {
						errMsg("algorithm", "addcut2Master", "ran out of memory to add new cut, include reduceCuts",0);
						return -1;
					}
                    dropIncumb = TRUE;
				}
                else {
					/* update the incumbent cut id */
					cell[i]->iCutIdx = cell[i]->cuts->cnt - 1;
                }

				/* add incumbent cut to the master problem in solver */
				cell[i]->cuts->val[cell[i]->iCutIdx]->alphaIncumb = cell[i]->cuts->val[cell[i]->iCutIdx]->alpha -
						vXv(cell[i]->cuts->val[cell[i]->iCutIdx]->beta, cell[0]->incumbX, NULL, lenX);
				status = addRow(cell[0]->sp->lp, lenX + 1, cell[i]->cuts->val[cell[i]->iCutIdx]->alphaIncumb, GE, 0, indices,
						cell[i]->cuts->val[cell[i]->iCutIdx]->beta);
				if (status){
					errMsg("solver", "addcut2Master", "failed to add new row to problem in solver", 0);
					return 1;
				}

				/* note the row number of the incumbent cut */
				cell[i]->cuts->val[cell[i]->iCutIdx]->rowNum = cell[0]->sp->mar++;

#ifdef CUT_CHECK
				writeProblem(cell[0]->sp->lp,"master_wIncumbCut.lp");
#endif
			}

            //MARK:check
            //cell[i]->pushcnt++;

			/* check to see if there is room for the candidate cut, else drop a cut */
			if (cell[i]->cuts->cnt > cell[i]->maxCuts) {
//            if (cell[i]->pushcnt > cell[i]->maxCuts) {
				/* status is the position of the oldest cut (which is dropped) and now holds the candidate cut */
				status = reduceCuts(cell, cell[0]->candidX, cell[0]->pi, cell[0]->lbType, lenX, lb, i, cell[0]->sp->lp);
				if ( status < 0 ) {
					errMsg("algorithm", "addCuts2Master", "failed to reduce cuts to make room for candidate cut", 0);
					return 1;
				}
				cut = cell[i]->cuts->val[status];
			}
			else {
                if (cell[i]->iCutUpdt == cell[0]->k) {
                    if (!(dropIncumb))
                        cut = cell[i]->cuts->val[cell[i]->cuts->cnt - 2];
                    else
                        /* when both candidate and incumbent cuts are added, it is the pen-ultimate cut */
                        cut = cell[i]->cuts->val[cell[i]->cuts->cnt - 1];
                }
				else
					/* when only candidate cut is added, it is the last cut */
					cut = cell[i]->cuts->val[cell[i]->cuts->cnt - 1];
			}

			/* add the row for candidate cut in the solver */
			cut->alphaIncumb = cut->alpha - vXv(cut->beta, cell[0]->incumbX, NULL, lenX);
			status = addRow(cell[0]->sp->lp, lenX + 1, cut->alphaIncumb, GE, 0, indices, cut->beta);
			if (status){
				errMsg("solver", "addcut2Master", "failed to add new row to problem in solver", 0);
				return 1;
			}

			/* note the row number of the candidate cut */
			cell[i]->cuts->val[cell[i]->cCutIdx]->rowNum = cell[0]->sp->mar++;
#ifdef CUT_CHECK
			writeProblem(cell[0]->sp->lp,"master_wCandidCut.lp");
#endif
		}
	}//END multi-cut version
	else {
        /* a single aggregated cut is added to the master problem */
        /* there is a single eta column in the master problem */
        indices[0] = lenX;

        /* compute the aggregated candidate cut */
        /* The istar field in aggregated cut is used to point to individual agent cuts */
        cut = newCut(lenX, numAgents, cell[0]->k);
        cut->cutObs = cell[0]->k;
        cut->isIncumb = FALSE;
        cut->beta[0] = 1.0;

        for ( i = 1; i < numAgents; i++ ) {
            cut->alpha += cell[i]->weight*cell[i]->cuts->val[cell[i]->cCutIdx]->alpha;
            for (c = 1; c <= lenX; c++)
                cut->beta[c] += cell[i]->weight*cell[i]->cuts->val[cell[i]->cCutIdx]->beta[c];
            cut->iStar[i] = cell[i]->cCutIdx;
        }
        cut->alphaIncumb = cut->alpha - vXv(cut->beta, cell[0]->incumbX, NULL, lenX);
        cell[0]->cuts->val[cell[0]->cuts->cnt] = cut;
        //cell[0]->cCutIdx = cell[0]->cuts->cnt;
        cell[0]->cuts->cnt++;

#ifdef CUT_CHECK
        printf("\tAggregated candidate cut: alpha = %lf, beta = ",  cut->alpha);
        printVector(cut->beta, lenX, NULL);
#endif

        /* check to see if there is room for the candidate cut, else drop a cut */
        if (cell[0]->cuts->cnt > cell[0]->maxCuts) {
            //            cell[0]->cuts->val[cell[0]->cuts->cnt-1]->rowNum = -1;

            i = 0;

            /* status is the position of the oldest cut (which is dropped) and now holds the candidate cut */
            status = reduceCuts(cell, cell[0]->candidX, cell[0]->pi, cell[0]->lbType, lenX, lb, i, cell[0]->sp->lp);
            if ( status < 0 ) {
                errMsg("algorithm", "addCuts2Master", "failed to add reduce cuts to make room for candidate cut", 0);
                return 1;
            }
            cell[0]->cCutIdx = status;
        }
        else {
            cell[0]->cCutIdx = cell[0]->cuts->cnt-1;
        }

        /* add the row for candidate cut to the master problem in solver */
        status = addRow(cell[0]->sp->lp, lenX + 1, cut->alphaIncumb, GE, 0, indices, cut->beta);
        if (status){
            errMsg("solver", "addcut2Master", "failed to add new row to problem in solver", 0);
            return 1;
        }

        /* note the row number of the candidate cut */
        cell[0]->cuts->val[cell[0]->cCutIdx]->rowNum = cell[0]->sp->mar++;
#ifdef CUT_CHECK
        writeProblem(cell[0]->sp->lp,"master_wCandidCut.lp");
#endif

        /* if an incumbent cut was created in the current iteration, then add the aggregated incumbent cut. */
        //MARK: check whether cell[0]->iCutUpdt updates
		if ((cell[0]->k - cell[0]->iCutUpdt) % config.TAU == 0) {
			/* compute the aggregated incumbent cut */
			/* The istar field in aggregated cut is used to point to individual agent cuts */
			cut = newCut(lenX, numAgents, cell[0]->k);
			cut->cutObs = cell[0]->k;	cut->isIncumb = TRUE;
			cut->beta[0] = 1.0;

			for ( i = 1; i < numAgents; i++ ) {
				cut->alpha = cell[i]->weight*cell[i]->cuts->val[cell[i]->iCutIdx]->alpha;
				for (c = 1; c <= lenX; c++)
					cut->beta[c] = cell[i]->weight*cell[i]->cuts->val[cell[i]->iCutIdx]->beta[c];
				cut->iStar[i] = cell[i]->iCutIdx;
			}
			cut->alphaIncumb = cut->alpha - vXv(cut->beta, cell[0]->incumbX, NULL, lenX);
			cell[0]->cuts->val[cell[0]->cuts->cnt] = cut;
			cell[0]->cuts->cnt++;

#ifdef CUT_CHECK
    printf("\t\t\t Aggregated incumbent cut: alpha = %lf, beta = ",  cut->alpha);
    printVector(cut->beta, lenX, NULL);
#endif

			if (cell[0]->cuts->cnt > cell[0]->maxCuts) {
//                cell[0]->cuts->val[cell[0]->cuts->cnt-1]->rowNum = -1;
                //i = 0;
				status = dropCut(cell, cell[0]->iCutIdx, 0, cell[0]->sp->lp);
				if ( status ) {
					errMsg("algorithm", "addcut2Master", "ran out of memory to add new cut, include reduceCuts",0);
					return -1;
				}
			}
			else
				/* update the incumbent cut id */
				cell[0]->iCutIdx = cell[0]->cuts->cnt - 1;

			status = addRow(cell[0]->sp->lp, lenX + 1, cut->alphaIncumb, GE, 0, indices, cut->beta);
			if (status){
				errMsg("solver", "addcut2Master", "failed to add new row to problem in solver", 0);
				return 1;
			}

			/* note the row number of the incumbent cut */
			cell[0]->cuts->val[cell[0]->iCutIdx]->rowNum = cell[0]->sp->mar++;

#ifdef CUT_CHECK
			writeProblem(cell[0]->sp->lp,"master_wIncumbCut.lp");
#endif
		}

	}

	mem_free(indices);
	return 0;
}//END addCuts2Master()

/* This function performs the updates on all the coefficients of eta in the master problem constraint matrix.  During every iteration,
 * each of the coefficients on eta are increased, so that the effect of the cut on the objective function is decreased. */
int changeEtaCol(LPptr lp, int numRows, int numCols, int k, cutType *cuts, double lb, int agentIdx, double weight) {
	double	etaCoef[1], etaBds[1], coef[1];
	int 	status, c, etaCol[1];
	char	bdsType[1];

	etaCol[0] = numCols;
	bdsType[0] = 'L';

	for (c = 0; c < cuts->cnt; c++){
		/* Currently both incumbent and candidate cuts are treated similarly, and sunk as iterations proceed */
		coef[0] = (double) (k) / (double) cuts->val[c]->cutObs;         // coeff k/j of eta column

		status = changeCol(lp, numCols + agentIdx -1, coef, cuts->val[c]->rowNum, cuts->val[c]->rowNum+1);
		if ( status ) {
			errMsg("solver", "chgEtaCol", "failed to change eta column in the stage problem", 0);
			return 1;
		}
	}

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

	return 0;
}//END chgEtaCol()

int updateRHS(cellType *master, int IniRow, int k) {
	int 	cnt;
	vector	rhs;
	intvec	indices;

	if (!(rhs = arr_alloc(master->sp->mar - IniRow, double)))
		errMsg("allocation", "updateRHS", "rhs", 0);

	if (!(indices = arr_alloc(master->sp->mar - IniRow, int)))
		errMsg("allocation", "updateRHS", "indices", 0);

	/* Now we change the right-hand of the master problem. */
	cnt = changeRHS(master->sp->lp, master->sp->mar - IniRow, indices, rhs);
	if (cnt)	{
		errMsg("solver", "changeQPrhs", "failed to change the right-hand side in the solver", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(indices);

	return 0;
}//END updateRHS

/* Construct the Q diagonal matrix and copy it for quadratic problem. */
int constructQP(LPptr lp, int numCols, double sigma) {
	int status = 0, idx, i;
	double *qsepvec;

	if (!(qsepvec = arr_alloc(numCols + numAgents - 1, double)))
		errMsg("Allocation", "constructQP", "qsepvec",0);

	/* Construct Q matrix, which is simply a diagonal matrix. */
	for (idx = 0; idx < numCols; idx++)
		qsepvec[idx] = 0.5 * sigma;

	/* This is for eta column */
	for (i = 0; i < numAgents - 1; i++)
		qsepvec[numCols + i] = 0.0;

	/* Now copy the Q matrix for QP problem. */
	status = copyQPseparable(lp, qsepvec);
	if (status) {
		fprintf(stderr, "Failed to copy Q matrix.\n");
		return 1;
	}

	mem_free(qsepvec);
	return 0;
}//END constructQP

int changeQPrhs(probType *prob, cellType **cell) {
	int 	cnt, n, status, offset;
	vector 	rhs;
	intvec 	indices;

    if (config.MULTI_CUT) {
        if (!(rhs =(vector) arr_alloc(cell[0]->sp->mar+numAgents*cell[0]->maxCuts+1, double)))
            errMsg("Allocation", "changeRhs", "rhs",0);
        if (!(indices =(intvec) arr_alloc(cell[0]->sp->mar+numAgents*cell[0]->maxCuts, int)))
            errMsg("Allocation", "changeRhs", "indices",0);

        /* Be careful with the one_norm!! In the CxX() routine, it assumes the 0th element is reserved for the 1_norm, in the returned vector, the T sparse
         vector, and the x vector. */
        for (cnt = 0; cnt < prob->num->rows; cnt++) {
            rhs[cnt + 1] = prob->sp->rhsx[cnt];
            indices[cnt] = cnt;
        }
        offset = prob->num->rows;

        /* b - A * xbar */
        rhs = MSparsexvSub(prob->Dbar, cell[0]->incumbX, rhs);

        /* change the right-hand side of individual cuts for every agent */
        for ( n = 1; n < numAgents; n++ ) {
            for ( cnt = 0; cnt < cell[n]->cuts->cnt; cnt++ ) {
                rhs[offset + cnt + 1] = cell[n]->cuts->val[cnt]->alpha - vXv(cell[n]->cuts->val[cnt]->beta, cell[0]->incumbX, NULL, prob->sp->mac);
                indices[offset + cnt] = cell[n]->cuts->val[cnt]->rowNum;

                cell[n]->cuts->val[cnt]->alphaIncumb = rhs[offset + cnt + 1];
            }
            offset += cell[n]->cuts->cnt;
        }
    }//END multi-cut version
    else {
        if (!(rhs =(vector) arr_alloc(cell[0]->sp->mar+cell[0]->maxCuts+1, double)))
            errMsg("Allocation", "changeRhs", "rhs",0);
        if (!(indices =(intvec) arr_alloc(cell[0]->sp->mar+cell[0]->maxCuts, int)))
            errMsg("Allocation", "changeRhs", "indices",0);

        /* Be careful with the one_norm!! In the CxX() routine, it assumes the 0th element is reserved for the 1_norm, in the returned vector, the T sparse
         vector, and the x vector. */
        for (cnt = 0; cnt < prob->num->rows; cnt++) {
            rhs[cnt + 1] = prob->sp->rhsx[cnt];
            indices[cnt] = cnt;
        }
        offset = prob->num->rows;

        /* b - A * xbar */
        rhs = MSparsexvSub(prob->Dbar, cell[0]->incumbX, rhs);

        /* change the right-hand side of individual cuts for every agent */
        for ( cnt = 0; cnt < cell[0]->cuts->cnt; cnt++ ) {
            rhs[offset + cnt + 1] = cell[0]->cuts->val[cnt]->alpha - vXv(cell[0]->cuts->val[cnt]->beta, cell[0]->incumbX, NULL, prob->sp->mac);
            indices[offset + cnt] = cell[0]->cuts->val[cnt]->rowNum;

            cell[0]->cuts->val[cnt]->alphaIncumb = rhs[offset + cnt + 1];
        }
        offset += cell[0]->cuts->cnt;
    }

#ifdef RHS_CHECK
	writeProblem(cell[0]->sp->lp, "QPchangeRHSBef.lp");
#endif
	/* Now we change the right-hand of the master problem. */
	status = changeRHS(cell[0]->sp->lp, offset, indices, rhs + 1);
	if ( status ) {
		errMsg("algorithm", "changeRhs", "failed to change the rhs", 0);
		return 1;
	}
#ifdef RHS_CHECK
	writeProblem(cell[0]->sp->lp, "QPchangeRHSAft.lp");
#endif

	mem_free(rhs);
	mem_free(indices);
	return 0;
}//END changeQPrh

/* This function changes the (lower) bounds of the variables, while changing from x to d. The lower bounds of d varibles are -xbar
 * (incumbent solution). */
int changeQPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector xk) {
	int 	status = 0, cnt;
	vector	lbounds, ubounds;
	intvec	lindices, uindices;
	char 	*llu, *ulu;

	if (!(lbounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "changeBounds", "lbounds",0);
	if (!(lindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "change_bounds", "lindices",0);
	if (!(llu = arr_alloc(numCols, char)))
		errMsg("Allocation", "changeBounds", "llu",0);

	if (!(ubounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "change_bounds", "ubounds",0);
	if (!(uindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "changeBounds", "uindices",0);
	if (!(ulu = arr_alloc(numCols, char)))
		errMsg("Allocation", "changeBounds", "ulu",0);

	/* Change the Upper Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		ubounds[cnt] = bdu[cnt] - xk[cnt + 1];
		uindices[cnt] = cnt;
		ulu[cnt] = 'U';
	}

	status = changeBDS(lp, numCols, uindices, ulu, ubounds);
	if (status) {
		errMsg("algorithm", "changeQP", "failed to change the upper bound in the solver", 0);
		return 1;
	}

	/* Change the Lower Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		lbounds[cnt] = bdl[cnt] - xk[cnt + 1];
		lindices[cnt] = cnt;
		llu[cnt] = 'L';
	}

	status = changeBDS(lp, numCols, lindices, llu, lbounds);
	if (status) {
		errMsg("algorithm", "changeQP", "failed to change the lower bound in the solver", 0);
		return 1;
	}

	mem_free(lbounds); mem_free(lindices); mem_free(llu);
	mem_free(ubounds); mem_free(uindices); mem_free(ulu);

	return 0;
}//END changeQPbds()

/**********************************************************************\
 ** This subroutine initializes the master problem by
 ** copying information from the decomposed prob[0](type: oneProblem)
 ** and adding a column for theta for modified benders decompostion.
     \**********************************************************************/
cellType *newMaster(probType *prob, vector xk, vector weight, double AggWeight) {
	cellType    *mcell;
    runTimeType *runTime;
	vector      rhs;
	intvec      indices;
	int         r, i, j, idx, cnt,status;
	long        colOffset, rowOffset;
	char        *q, etaName[NAMESIZE];

	if (!(mcell = (cellType *) mem_malloc (sizeof(cellType))))
		errMsg("Memory allocation", "new_master", "Faile to allocate memory to master cell", 0);
	mcell->ID = 0;
	mcell->weight = AggWeight;

	if (!(mcell->sp = (oneProblem *) mem_malloc (sizeof(oneProblem))))
		errMsg("Memory allocation", "new_master", "Faile to allocate memory to mcell->sp", 0);

	/* -+-+-+-+-+-+-+-+-+-+-+-+-+-+- Allocating memory to mcell-sp -+-+-+-+-+-+-+-+-+-+-+-+-+-+- */
	mcell->sp->type = config.MASTERTYPE;                  /* type of problem: LP, QP, MIP or MIQP */
	mcell->sp->objsen = prob->sp->objsen;                 /* sense of the objective: 1 for minimization and -1 for maximization */
	mcell->sp->mar = prob->sp->mar;                       /* number of rows */
	mcell->sp->numInt = prob->sp->numInt;                 /* number of integer variables in the problem  */
	mcell->sp->numnz = prob->sp->numnz;                   /* number of non-zero elements in constraint matrix */
	mcell->sp->matsz = prob->sp->matsz;                   /* extended matrix size */
	mcell->sp->marsz = prob->sp->marsz;                   /* extended row size */
	mcell->sp->rstorsz = prob->sp->rstorsz;               /* memory size for storing row names */
	if ( config.MULTI_CUT ) {
		mcell->sp->mac 		= prob->sp->mac+numAgents-1;           /* number of columns + etas */
		mcell->sp->macsz 	= prob->sp->macsz + numAgents - 1;       /* extended column size */
		mcell->sp->cstorsz 	= prob->sp->cstorsz + (numAgents - 1) * NAMESIZE;    /* memory size for storing column names */
	}
	else {
		mcell->sp->mac 		= prob->sp->mac+1;           		/* number of columns + etas */
		mcell->sp->macsz 	= prob->sp->macsz + 1;       		/* extended column size */
		mcell->sp->cstorsz 	= prob->sp->cstorsz + NAMESIZE;    	/* memory size for storing column names */
	}

	/* Allocate memory to the information whose type is string */
	if (!(mcell->sp->name = (string) arr_alloc(NAMESIZE, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to mcell->sp->name",0);
	if (!(mcell->sp->senx = (string) arr_alloc(mcell->sp->marsz,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to mcell->sp->senx",0);
	if (!(mcell->sp->ctype = (string) arr_alloc(mcell->sp->macsz,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to mcell->sp->ctype",0);
	if (!(mcell->sp->objname = (string) arr_alloc(NAMESIZE,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to mcell->sp->objname",0);
	if (!(mcell->sp->rname = (string *) arr_alloc(mcell->sp->marsz,string)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to mcell->sp->rname",0);
	if (!(mcell->sp->rstore = (string) arr_alloc(mcell->sp->rstorsz, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to mcell->sp->rstore",0);
	if (!(mcell->sp->cname = (string*) arr_alloc(mcell->sp->macsz,string)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to mcell->sp->cname",0);
	if (!(mcell->sp->cstore = (string) arr_alloc(mcell->sp->cstorsz, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to mcell->sp->cstore",0);

	/* Allocate memory to the information whose type is vector */
	if (!(mcell->sp->objx = (vector) arr_alloc(mcell->sp->macsz, double)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to mcell->sp->objx",0);
	if (!(mcell->sp->rhsx = (vector) arr_alloc(mcell->sp->marsz, double)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to mcell->sp->rhsx",0);
	if (!(mcell->sp->matval = (vector) arr_alloc(mcell->sp->matsz, double)))
		errMsg("allocation", "new_master", "mcell->sp->matval",0);
	if (!(mcell->sp->bdl = (vector) arr_alloc(mcell->sp->macsz, double)))
		errMsg("allocation", "new_master", "mcell->sp->bdl",0);
	if (!(mcell->sp->bdu = (vector) arr_alloc(mcell->sp->macsz, double)))
		errMsg("allocation", "new_master", "mcell->sp->bdu",0);

	/* Allocate memory to the information whose type is intvec */
	if (!(mcell->sp->matbeg = (intvec) arr_alloc(mcell->sp->macsz, int)))
		errMsg("allocation", "new_master", "mcell->sp->matbeg",0);
	if (!(mcell->sp->matcnt = (intvec) arr_alloc(mcell->sp->macsz, int)))
		errMsg("allocation", "new_master", "mcell->sp->matcnt",0);
	if (!(mcell->sp->matind = (intvec) arr_alloc(mcell->sp->matsz, int)))
		errMsg("allocation", "new_master", "mcell->sp->matind",0);

	strcpy(mcell->sp->name, prob->sp->name);           /* Copy problem name */
	strcpy(mcell->sp->objname, prob->sp->objname);     /* Copy objective name */

	/* Copy problem's column and row names */
	i = 0;
	for (q = prob->sp->cname[0]; q < prob->sp->cname[0] + prob->sp->cstorsz; q++)
		mcell->sp->cstore[i++] = *q;

	i = 0;
	for (q = prob->sp->rname[0]; q < prob->sp->rname[0] + prob->sp->rstorsz; q++)
		mcell->sp->rstore[i++] = *q;

	/* Calculate difference in pointers for master/copy row and column names */
	colOffset = mcell->sp->cstore - prob->sp->cname[0];
	rowOffset = mcell->sp->rstore - prob->sp->rname[0];

	/* Copy the all column information from the original master problem */
	cnt = 0;
	for (j = 0; j < prob->sp->mac; j++) {
		/* Copy objective function coefficients */
		mcell->sp->objx[j] = prob->sp->objx[j];
		/* Copy the decision variable type */
		mcell->sp->ctype[j] = prob->sp->ctype[j];
		/* Copy the upper bound and lower bound */
		mcell->sp->bdu[j] = prob->sp->bdu[j];
		mcell->sp->bdl[j] = prob->sp->bdl[j];
		/* Copy column names, offset by length */
		mcell->sp->cname[j] = prob->sp->cname[j] + colOffset;
		/* Copy the master sparse matrix beginning position of each column */
		mcell->sp->matbeg[j] = cnt;
		/* Copy the sparse matrix non-zero element count */
		mcell->sp->matcnt[j] = prob->sp->matcnt[j];
		mcell->sp->ctype[j] = prob->sp->ctype[j];
		/* Loop through all non-zero elements in this column */
		for (idx = prob->sp->matbeg[j]; idx < prob->sp->matbeg[j] + prob->sp->matcnt[j]; idx++) {
			/* Copy the non-zero coefficient */
			mcell->sp->matval[cnt] = prob->sp->matval[idx];
			/* Copy the row entry of the non-zero elements */
			mcell->sp->matind[cnt] = prob->sp->matind[idx];
			cnt++;
		}
	}

	/* Copy all information concerning rows of master */
	for (r = 0; r < prob->sp->mar; r++) {
		/* Copy the right hand side value */
		mcell->sp->rhsx[r] = prob->sp->rhsx[r];
		/* Copy the constraint sense */
		mcell->sp->senx[r] = prob->sp->senx[r];
		/* Copy row names, offset by length */
		mcell->sp->rname[r] = prob->sp->rname[r] + rowOffset;
	}

	/* Initialize information for the extra column in the new master. */
	colOffset = prob->sp->cstorsz;
	if ( config.MULTI_CUT )
		for (i = 1; i < numAgents; i++) {
			sprintf(etaName, "eta%d",i);
			strcat(mcell->sp->cstore + colOffset, etaName);
			mcell->sp->cname[prob->sp->mac - 1 + i] = mcell->sp->cstore + colOffset;
			mcell->sp->objx[prob->sp->mac - 1 + i] = weight[i];			// prob->sp->mac is the last column in the original master
			mcell->sp->ctype[prob->sp->mac - 1 + i] = 'C';
			mcell->sp->bdu[prob->sp->mac - 1 + i] = INFBOUND;
			mcell->sp->bdl[prob->sp->mac - 1 + i] = 0.0;
			mcell->sp->matbeg[prob->sp->mac - 1 + i] = prob->sp->numnz;	// Beginning point in matval/matind in eta columns. every eta column begins at the same address
			mcell->sp->matcnt[prob->sp->mac - 1 + i] = 0;               // Only optimality cuts has eta
			colOffset += strlen(etaName)+1;
		}
	else {
		strcpy(mcell->sp->cstore + prob->sp->cstorsz, "eta");
		mcell->sp->cname[prob->sp->mac] = mcell->sp->cstore + colOffset;
		mcell->sp->objx[prob->sp->mac] = 1.0;			// prob->sp->mac is the last column in the original master
		mcell->sp->ctype[prob->sp->mac] = 'C';
		mcell->sp->bdu[prob->sp->mac] = INFBOUND;
		mcell->sp->bdl[prob->sp->mac] = 0.0;
		mcell->sp->matbeg[prob->sp->mac] = prob->sp->numnz;	// Beginning point in matval/matind in eta columns. every eta column begins at the same address
		mcell->sp->matcnt[prob->sp->mac] = 0;               // Only optimality cuts has eta
	}

	/* Load the copy into CPLEX */
	mcell->sp->lp = setupProblem(mcell->sp->name, mcell->sp->type, mcell->sp->mac, mcell->sp->mar, mcell->sp->objsen, mcell->sp->objx, mcell->sp->rhsx, mcell->sp->senx, mcell->sp->matbeg, mcell->sp->matcnt,mcell->sp->matind, mcell->sp->matval, mcell->sp->bdl, mcell->sp->bdu, NULL, mcell->sp->cname, mcell->sp->rname, mcell->sp->ctype);
	if ( mcell->sp->lp == NULL ) {
		errMsg("Problem Setup", "new_master", "failed to setup master problem in the solver",0);
		return NULL;
	}

#if 0
	status = writeProblem(mcell->sp->lp, "newMaster.lp");
	if ( status ) {
		errMsg("write problem", "new_master", "failed to write master problem to file",0);
		return NULL;
	}
#endif

	/* -+-+-+-+-+-+-+-+-+-+-+ Allocating memory to other variables that belongs to master mcell +-+-+-+-+-+-+-+-+-+- */
	//MARK: original onfig.CUT_MULT * MAX_CUTS(prob->num->cols),  MAX_CUTS = (3 * (c) + 3)
	if (prob->lb == 0)
		mcell->lbType = TRIVIAL;
	else
		mcell->lbType = NONTRIVIAL;

	if ( config.MULTI_CUT ) {
		mcell->maxCuts = config.CUT_MULT * MAX_CUTS(prob->num->cols);
		mcell->cuts = NULL;
	}
	else {
		mcell->maxCuts = config.CUT_MULT * MAX_CUTS(prob->num->cols);
		mcell->cuts = newCuts(mcell->maxCuts);
	}
	mcell->candidX = duplicVector(xk, prob->num->cols);
	mcell->candidEst = prob->lb + vXvSparse(mcell->candidX,prob->dBar);
	mcell->pi_ratio = NULL;
	mcell->dualStableFlag = FALSE;
	mcell->fcuts = NULL;
	mcell->incumbEst = mcell->candidEst;
	mcell->quadScalar = config.MIN_QUAD_SCALAR;     /* The quadratic scalar, 'sigma'*/
	mcell->LPcnt = 0;
	mcell->feasCnt = 0;
	mcell->infeasIncumb = FALSE;
	mcell->feasFlag = TRUE;
	mcell->incumbStdev = 0.0;
	mcell->incumbChg = FALSE;
	mcell->iCutIdx = 0;
	mcell->iCutUpdt = 0;
	mcell->gamma = 0.0;
	mcell->normDk_1 = 0.0;
	mcell->normDk = 0.0;
	mcell->optValM = 0.0;
	mcell->optFlag = FALSE;
	mcell->k = 0;
	mcell->spRHS = NULL;
	mcell->full_test_error = 0.0;

	/* solution parts of the cell */
	if ( !(mcell->pi = (vector) arr_alloc(prob->num->rows + numAgents*(mcell->maxCuts + 1), double)) )
		errMsg("allocation", "newMaster", "mcell->pi", 0);
	if ( !(mcell->di = (vector) arr_alloc(prob->num->cols + 2, double)) )
		errMsg("allocation", "newMaster", "mcell->di", 0);

	if (!(rhs =(vector) arr_alloc(prob->num->rows+ mcell->maxCuts + 1, double)))
		errMsg("Allocation", "newMaster", "rhs",0);
	if (!(indices =(intvec) arr_alloc(prob->num->rows+ mcell->maxCuts + 1, int)))
		errMsg("Allocation", "newMaster", "indices",0);

	if (config.MASTERTYPE == PROB_QP) {
		mcell->incumbX = duplicVector(xk, prob->num->cols);

		/* change QP rhs */
		/* Be careful with the one_norm!! In the CxX() routine, it assumes the 0th element is reserved for the 1_norm, in the returned vector, the T sparse
         vector, and the x vector. */
		for (cnt = 0; cnt < prob->num->rows; cnt++) {
			rhs[cnt + 1] = prob->sp->rhsx[cnt];
			indices[cnt] = cnt;
		}

		/* b - A * xbar */
		rhs = MSparsexvSub(prob->Dbar, mcell->incumbX, rhs);

#ifdef RHS_CHECK
		writeProblem(mcell->sp->lp, "NewQPchangeRHSBef.lp");
#endif
		/* Now we change the right-hand of the master problem. */
		status = changeRHS(mcell->sp->lp, prob->num->rows, indices, rhs + 1);
		if ( status ) {
			errMsg("algorithm", "newMaster", "failed to change the rhs", 0);
			return NULL;
		}
#ifdef RHS_CHECK
		writeProblem(mcell->sp->lp, "NewQPchangeRHSAft.lp");
#endif

		/* change QP bounds */
		status = changeQPbds(mcell->sp->lp, prob->num->cols, prob->sp->bdl, prob->sp->bdu, mcell->incumbX);
		if ( status ) {
			errMsg("algorithm", "newMaster", "failed to change the bounds", 0);
			return NULL;
		}
	}
	else
		mcell->incumbX = NULL;

	mcell->sigma = NULL;
	mcell->delta = NULL;
	mcell->omega = NULL;
	mcell->lambda = NULL;

    /* initialize the time structure */
    if (!(mcell->runTime =(runTimeType *) mem_malloc(sizeof(runTimeType))))
        errMsg("Allocation", "newMaster", "runTime", 0);
    mcell->runTime->iterCutGenTime = 0.0; mcell->runTime->iterSolTime = 0.0; mcell->runTime->iterTime = 0.0;
    mcell->runTime->totCutGenTime = 0.0; mcell->runTime->totSolTime = 0.0; mcell->runTime->totTime = 0.0;

	mem_free(rhs);
	mem_free(indices);
	return mcell;
}//END newMaster
