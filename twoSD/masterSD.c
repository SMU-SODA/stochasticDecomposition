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

#include "twoSD.h"

extern configType config;

/* This function is the regularized QP version of master problem. The master problem is solved after the newest cut is added to master problem,
 the incumbent cut is updated if necessary. Here the coefficients on all the cuts are updated, and finally master problem is solved. */
int solveSDMaster(numType *num, sparseVector *dBar, cellType *cell, int IniRow, double lb) {
	double 	d2 = 0.0; /* height at the candidate solution. */
	int 	status, i;

	if( changeEtaCol(cell->master->lp, num->rows, num->cols, cell->k, cell->cuts, lb) ) {
		errMsg("algorithm", "solveMaster", "failed to change the eta column coefficients", 0);
		return 1;
	}

	if ( cell->lbType == NONTRIVIAL ) {
		/* update the right-hand side of cuts to reflect the non-trivial lower bound */
		if ( updateRHS(cell->master->lp, cell->cuts, cell->k, cell->lb) ) {
			errMsg("algorithm", "solveQPMaster", "failed to update right-hand side with lower bound information", 0);
			return 1;
		}
	}

#ifdef ALGO_CHECK
	writeProblem(cell->master->lp, "masterCell.lp");
#endif

	/* solve the master problem */
	if ( solveProblem(cell->master->lp, cell->master->name, config.MASTER_TYPE, &status) ) {
		writeProblem(cell->master->lp, "error.lp");
		errMsg("algorithm", "solveMaster", "failed to solve the master problem", 0);
		return 1;
	}

	/* increment the number of problems solved during algorithm */
	cell->LPcnt++;

	/* Get the most recent optimal solution to master program */
	if ( getPrimal(cell->master->lp, cell->candidX, num->cols) ) {
		errMsg("algorithm", "solveMaster", "failed to obtain the primal solution for master", 0);
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

	/* Get the dual solution too */
	if ( getDual(cell->master->lp, cell->piM, cell->master->mar) ) {
		errMsg("solver", "solveQPMaster", "failed to obtain dual solutions to master", 0);
		return 1;
	}
	if ( getDualSlacks(cell->master->lp, cell->djM, num->cols) ) {
		errMsg("solver", "solveQPMaster", "failed to obtain dual slacks for master", 0);
		return 1;
	}

	/* Find the highest cut at the candidate solution. where cut_height = alpha - beta(xbar + \Delta X) */
	cell->candidEst = vXvSparse(cell->candidX, dBar) + maxCutHeight(cell->cuts, cell->candidX, num->cols, TRUE, cell->k, lb);

	/* Calculate gamma for next improvement check on incumbent x. */
	cell->gamma = cell->candidEst - cell->incumbEst;

	return 0;
}//END solveSDMaster()

/* This function performs the updates on all the coefficients of eta in the master problem constraint matrix.  During every iteration,
 * each of the coefficients on eta are increased, so that the effect of the cut on the objective function is decreased. */
int changeEtaCol(LPptr lp, int numRows, int numCols, int k, cutsType *cuts, double lb) {
	double	etaCoef[1], etaBds[1], coef[1];
	int 	status, c, etaCol[1];
	char	bdsType[1];

	etaCol[0] = numCols;
	bdsType[0] = 'L';

	for (c = 0; c < cuts->cnt; c++){
		/* Currently both incumbent and candidate cuts are treated similarly, and sunk as iterations proceed */
		coef[0] = (double) (k) / (double) cuts->vals[c]->numSamples;         // coeff k/j of eta column

		status = changeCol(lp, numCols, coef, cuts->vals[c]->rowNum, cuts->vals[c]->rowNum+1);
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

int updateRHS(LPptr lp, cutsType *cuts, int numIter, double lb) {
	int 	cnt;
	vector	rhs;
	intvec	indices;

	if (!(rhs = arr_alloc(cuts->cnt, double)))
		errMsg("allocation", "updateRHS", "rhs", 0);
	if (!(indices = arr_alloc(cuts->cnt, int)))
		errMsg("allocation", "updateRHS", "indices", 0);

	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		rhs[cnt] = cuts->vals[cnt]->alphaIncumb + ((double) numIter / (double) cuts->vals[cnt]->numSamples - 1) * lb;
		indices[cnt] = cuts->vals[cnt]->rowNum;
	}

	/* Now we change the right-hand of the master problem. */
	if ( changeRHS(lp, cuts->cnt, indices, rhs) ) {
		errMsg("solver", "changeQPrhs", "failed to change the right-hand side in the solver", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(indices);

	return 0;
}//END updateRHS

int formSDCut(probType *prob, cellType *cell, vector Xvect, int omegaIdx, BOOL newOmegaFlag, BOOL isIncumb) {
	oneCut *cut;
	int    cutIdx;

	/* (a) Construct the subproblem with input observation and master solution, solve the subproblem, and complete stochastic updates */
	if ( solveSubprob(prob, cell->subprob, Xvect, cell->basis, cell->lambda, cell->sigma, cell->delta, config.MAX_ITER,
			cell->omega, omegaIdx, newOmegaFlag, cell->k, config.TOLERANCE) ) {
		errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
		return -1;
	}

	/* (b) create an affine lower bound */
	cut = SDCut(prob->num, prob->coord, cell->basis, cell->sigma, cell->delta, cell->omega, Xvect, cell->k, &cell->dualStableFlag, cell->pi_ratio, cell->lb);
	if ( cut == NULL ) {
		errMsg("algorithm", "formSDCut", "failed to create the affine minorant", 0);
		return -1;
	}

	/* (c) add cut to the master problem  */
	if ( (cutIdx = addCut2Master(cell, cut, TRUE, prob->num->prevCols, cell->lb)) < 0 ) {
		errMsg("algorithm", "formSDCut", "failed to add the new cut to master problem", 0);
		return -1;
	}
	return cutIdx;
}//END formCut()

oneCut *SDCut(numType *num, coordType *coord, basisType *basis, sigmaType *sigma, deltaType *delta, omegaType *omega, vector Xvect, int numSamples,
		BOOL *dualStableFlag, vector pi_ratio, double lb) {
	oneCut *cut;
	vector 	piCbarX, beta;
	double  argmaxAll, argmaxNew, argmax, alpha = 0.0, argmax_dif_sum = 0.0, argmax_all_sum = 0.0, variance = 1.0, mean;
	int 	istarAll, istarNew, istar, c, cnt, obs, offset;
	BOOL    pi_eval_flag = FALSE;

	/* allocate memory to hold a new cut */
	cut = newCut(num->prevCols, omega->cnt, numSamples);

	/* Need to store  Pi x Cbar x X independently of observation loop */
	if (!(piCbarX= arr_alloc(sigma->cnt, double)))
		errMsg("Allocation", "SDCut", "pi_Tbar_x",0);
	if ( !(beta = (vector) arr_alloc(num->prevCols + 1, double)) )
		errMsg("Allocation", "SDCut", "beta", 0);

	/* Calculate pi_eval_flag to determine the way of computing argmax */
	if (numSamples > config.PI_EVAL_START && !(numSamples % config.PI_CYCLE))
		pi_eval_flag = TRUE;

	/* Calculate (Pi x Cbar) x X by mult. each VxT by X, one at a time */
	for (cnt = 0; cnt < sigma->cnt; cnt++)
		piCbarX[cnt] = vXv(sigma->vals[cnt].piC, Xvect, coord->colsC, num->cntCcols);

	offset = num->rvbOmCnt + num->rvCOmCnt;
	/* Test for omega issues */
	for (obs = 0; obs < omega->cnt; obs++) {
		/* For each observation, find the Pi which maximizes height at X. */
		if (pi_eval_flag == TRUE) {
			istarAll = computeIstar(num, coord, basis, sigma, delta, Xvect, piCbarX, omega->vals[obs]+offset, obs, numSamples, pi_eval_flag, &argmaxAll, FALSE);
			istarNew = computeIstar(num, coord, basis, sigma, delta, Xvect, piCbarX, omega->vals[obs]+offset, obs, numSamples, TRUE, &argmaxNew, TRUE);

			if (argmaxNew > argmaxAll) {
				argmax = argmaxNew; istar = istarNew;
			}
			else {
				argmax = argmaxAll; istar = istarAll;
			}

			argmax_dif_sum += max(argmaxAll - lb, 0) * omega->weights[obs];
			argmax_all_sum += max(argmax - lb, 0) * omega->weights[obs];
		}
		else {
			/* identify the maximal Pi for each observation */
			istar = computeIstar(num, coord, basis, sigma, delta, Xvect, piCbarX, omega->vals[obs], obs, numSamples, pi_eval_flag, &argmax, FALSE);
		}

		if ( istar < 0 ) {
			errMsg("algorithm", "SDCut", "failed to identify maximal Pi for an observation", 0);
			return NULL;
		}
		cut->iStar[obs] = istar;

		/* Average using these Pi's to calculate the cut itself (update alpha and beta) */
		alpha += (sigma->vals[basis->vals[istar]->sigmaIdx[0]].pib + delta->vals[basis->vals[istar]->lambdaIdx[0]][obs].pib)* omega->weights[obs];

		for (c = 1; c <= num->cntCcols; c++)
			beta[coord->colsC[c]] += sigma->vals[basis->vals[istar]->sigmaIdx[0]].piC[c] * omega->weights[obs];
		for (c = 1; c <= num->rvCOmCnt; c++)
			beta[coord->rvCols[c]] += delta->vals[basis->vals[istar]->lambdaIdx[0]][obs].piC[c] * omega->weights[obs];

		for ( cnt = 1; cnt <= basis->vals[istar]->phiLength; cnt++ ) {
			alpha += (sigma->vals[basis->vals[istar]->sigmaIdx[cnt]].pib + delta->vals[basis->vals[istar]->lambdaIdx[cnt]][obs].pib)
					* (omega->vals[obs][offset+basis->vals[istar]->omegaIdx[cnt]]) * omega->weights[obs];
			for (c = 1; c <= num->cntCcols; c++)
				beta[coord->colsC[c]] += sigma->vals[basis->vals[istar]->sigmaIdx[cnt]].piC[c] * (omega->vals[obs][offset+basis->vals[istar]->omegaIdx[cnt]]) * omega->weights[obs];
			for (c = 1; c <= num->rvCOmCnt; c++)
				beta[coord->rvCols[c]] += delta->vals[basis->vals[istar]->lambdaIdx[cnt]][obs].piC[c] * (omega->vals[obs][offset+basis->vals[istar]->omegaIdx[cnt]]) * omega->weights[obs];
		}
	}

	if (pi_eval_flag == TRUE) {
		pi_ratio[numSamples % config.SCAN_LEN] = argmax_dif_sum / argmax_all_sum;
		if (numSamples - config.PI_EVAL_START > config.SCAN_LEN)
			calcMeanVariance(pi_ratio, config.SCAN_LEN, &mean, &variance);

		if (DBL_ABS(variance) >= .000002 || (pi_ratio[numSamples % config.SCAN_LEN]) < 0.95)
			*dualStableFlag = FALSE;
		else
			*dualStableFlag = TRUE;
	}

	cut->alpha = alpha / numSamples;

	for (c = 1; c <= num->prevCols; c++)
		cut->beta[c] = beta[c] / numSamples;

	/* coefficient of eta coloumn */
	cut->beta[0] = 1.0;

	mem_free(piCbarX);
	mem_free(beta);

	return cut;
}//END SDCut

/* This function determines whether the "stagewise descent property" is satisified.  If the current approximation of f_k gives
 * a lower difference between the candidate and incumbent x than the previous approximation gave, then the incumbent x is
 * updated to the candidate x, and the reference to the incumbent cut is updated as well.  The function returns TRUE if the
 * incumbent was updated; FALSE otherwise. */
int checkImprovement(probType *prob, cellType *cell, int candidCut) {
	double  candidEst;

	/* Calculate height at new candidate x with newest cut included */
	candidEst = vXvSparse(cell->candidX, prob->dBar) + maxCutHeight(cell->cuts, cell->candidX, prob->num->cols, TRUE, cell->k, cell->lb);
	cell->incumbEst = vXvSparse(cell->incumbX, prob->dBar) + maxCutHeight(cell->cuts, cell->incumbX, prob->num->cols, TRUE, cell->k, cell->lb);

#ifdef SOL_CHECK
	printf("AggcandidEst =%lf, AggIncumEst =%lf\n",AggcandidEst, cell->incumbEst);
#endif

	/* If we see considerable improvement, then change the incumbent */
	if ((candidEst - cell->incumbEst) < (config.R1 * cell->gamma)) {
		/* when we find an improvement, then we need to replace the incumbent x with candidate x */
		if ( replaceIncumbent(prob, cell, candidEst) ) {
			errMsg("algorithm", "checkImprovement", "failed to replace incumbent solution with candidate", 0);
			return 1;
		}
		cell->iCutIdx = candidCut;
		cell->incumbChg = FALSE;
		printf("+"); fflush(stdout);
	}
	else {
		/* Update quad_scalar when no incumbent is found. */
		cell->quadScalar = min(config.MAX_QUAD_SCALAR, cell->quadScalar / config.R2);
		cell->normDk_1 = cell->normDk;

		/* change the proximal term in the solver */
		if ( changeQPproximal(cell->master->lp, prob->num->cols, cell->quadScalar) ) {
			errMsg("setup", "newCell", "failed to add the proximal term to QP", 0);
			return 1;
		}
	}

	return 0;
}//END checkImprovement()

int replaceIncumbent(probType *prob, cellType *cell, double candidEst) {

	/* replace the incumbent solution with the candidate solution */
	copyVector(cell->candidX, cell->incumbX, prob->num->cols, 1);
	cell->incumbEst = candidEst;

	/* update the proximal parameter based on estimated improvement */
	if ( cell->normDk > config.TOLERANCE )
		if ( cell->normDk >= config.R3 * cell->normDk_1 ) {
			cell->quadScalar *= config.R2 * config.R3 * cell->normDk_1/ cell->normDk;
			cell->quadScalar  = min(config.MAX_QUAD_SCALAR, cell->quadScalar);
			cell->quadScalar = max(config.MIN_QUAD_SCALAR, cell->quadScalar);
		}

	/* update the right-hand side and the bounds with new incumbent solution */
	if ( constructQP(prob, cell, cell->incumbX, cell->quadScalar) ) {
		errMsg("algorithm", "replaceIncumbent", "failed to change the right-hand side after incumbent change", 0);
		return 1;
	}

	/* update the candidate cut as the new incumbent cut */
	cell->iCutUpdt = cell->k;
	cell->incumbChg = TRUE;

	/* keep the two norm of solution*/
	cell->normDk_1 = cell->normDk;
	/* Since incumbent solution is now replaced by a candidate, we assume it is feasible now */
	cell->infeasIncumb = FALSE;
	/* gamma needs to be reset to 0 since there's no difference between candidate and incumbent*/
	cell->gamma = 0.0;

	return 0;
}//END replaceIncumbent()
