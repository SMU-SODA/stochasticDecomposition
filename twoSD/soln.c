/*
 * soln.c
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

/***********************************************************************\
 ** This function determines whether the "stagewise descent property" is
 ** satisified.  If the current approximation of f_k gives a lower difference
 ** between the candidate and incumbent x than the previous approximation
 ** gave, then the incumbent x is updated to the candidate x, and the
 ** reference to the incumbent cut is updated as well.  The function returns
 ** TRUE if the incumbent was updated; FALSE otherwise.
 \***********************************************************************/
BOOL checkImprovement(probType *prob, cellType **cell) {
	intvec  maxCuts;
    double  AggcandidEst;
    int     i, temp;

    return FALSE;
}
//
//	if (!(maxCuts = (intvec) arr_alloc(numAgents, int)))
//		errMsg("Allocation", "checkImprovement", "Fail to allocate memory to maxCuts",0);
//
//	//cell[0]->candidEst = vXvSparse(cell[0]->candidX, prob->dBar);
//    AggcandidEst = vXvSparse(cell[0]->candidX, prob->dBar);
//	cell[0]->incumbEst = vXvSparse(cell[0]->incumbX, prob->dBar);
//
//
//	AggcandidEst += maxCutHeight(cell[0]->lbType, cell[0]->cuts, cell[0]->k, cell[0]->candidX, prob->num->cols, prob->lb, &maxCuts[0]);
//	cell[0]->incumbEst += maxCutHeight(cell[0]->lbType, cell[0]->cuts, cell[0]->k, cell[0]->incumbX, prob->num->cols, prob->lb, &maxCuts[0]);
//#ifdef SOL_CHECK
//    printf("AggcandidEst =%lf, AggIncumEst =%lf\n",AggcandidEst, cell[0]->incumbEst);
//#endif
//	/* If we see considerable improvement, then change the incumbent */
//	if ((AggcandidEst - cell[0]->incumbEst) < (config.R1 * cell->gamma)) {
//		/* when we find an improvement, then we need to replace the incumbent x with candidate x */
//		replaceIncumbent(prob, cell, cell[0]->k, maxCuts, AggcandidEst);
//
//		mem_free(maxCuts);
//		return TRUE;
//	}
//	else {
//		/* Update quad_scalar when no incumbent is found. */
//		cell[0]->quadScalar = min(config.MAX_QUAD_SCALAR, cell[0]->quadScalar / config.R2);
//		cell[0]->incumbStdev *= (cell[0]->k - 1) / (double) (cell[0]->k);
//		cell[0]->normDk_1 = cell[0]->normDk;
//
//		mem_free(maxCuts);
//		return FALSE;
//	}
//}//END checkImprovement()
//
//int replaceIncumbent(probType *prob, cellType **cell, int k, intvec maxCutID, double AggcandidEst) {
//	int     i, status;
//
//	/* replace the incumbent solution with the candidate solution */
//	copyVector(cell[0]->candidX, cell[0]->incumbX, prob->num->cols, 1);
//	cell[0]->incumbEst = AggcandidEst;
//
//	/* update the proximal parameter based on estimated improvement */
//	if ( k > 1 && cell[0]->normDk > config.TOLERANCE )
//		if ( cell[0]->normDk >= config.R3 * cell[0]->normDk_1 ) {
//			cell[0]->quadScalar *= config.R2 * config.R3 * cell[0]->normDk_1/ cell[0]->normDk;
//			cell[0]->quadScalar  = min(config.MAX_QUAD_SCALAR, cell[0]->quadScalar);
//			cell[0]->quadScalar = max(config.MIN_QUAD_SCALAR, cell[0]->quadScalar);
//		}
//
//	/* update the right-hand side and the bounds with new incumbent solution */
//	status = changeQPrhs(prob, cell);
//	if ( status ) {
//		errMsg("algorithm", "replaceIncumbent", "failed to change the right-hand side after incumbent change", 0);
//		return 1;
//	}
//	status = changeQPbds(cell->master->lp, prob->num->cols, prob->sp->bdl, prob->sp->bdu, cell[0]->incumbX);
//	if ( status ) {
//		errMsg("algorithm", "replaceIncumbent", "failed to change the bounds after incumbent update", 0);
//		return 1;
//	}
//
//	/* update the candidate cut as the new incumbent cut */
//	cell[0]->iCutUpdt = cell[0]->k;
//	cell[0]->iCutIdx = cell[0]->cCutIdx;
//
//	cell[0]->incumbChg = TRUE;
//
//	/* keep the two norm of solution*/
//	cell[0]->normDk_1 = cell[0]->normDk;
//	/* Since incumbent solution is now replaced by a candidate, we assume it is feasible now */
//	cell[0]->infeasIncumb = FALSE;
//	/* gamma needs to be reset to 0 since there's no difference between candidate and incumbent*/
//	cell[0]->gamma = 0.0;
//
//	printf("(+%d)", cell[0]->k); fflush(stdout);
//
//	return 0;
//}//END replaceIncumbent()

/* This function loops through a set of cuts and find the highest cut height at the specified position x */
double maxCutHeight(int lbType, cutsType *cuts, int currIter, vector xk, int betaLen, double lb, int *maxCutID) {
	double Sm = -INF, ht = 0.0;
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		ht = cutHeight(lbType, cuts->vals[cnt], currIter, xk, betaLen, lb);
		if (Sm < ht) {
			Sm = ht;
			(*maxCutID) = cnt;
		}
	}

	return Sm;
}//END maxCutHeight

/* This function calculates and returns the height of a given cut at a given X.  It includes the k/(k-1) update, but does not include
 * the coefficients due to the cell. */
double cutHeight(int lbType, oneCut *cut, int currIter, vector xk, int betaLen, double lb) {
	double height;
	double t_over_k = ((double) cut->cutObs / (double) currIter);

	/* A cut is calculated as alpha - beta x X */
	height = cut->alpha - vXv(cut->beta, xk, NULL, betaLen);

	/* Weight cut based on number of observations used to form it */
	height *= t_over_k;

	/* Updated for optimality cut height*/
	if (lbType == NONTRIVIAL)
		height += (1 - t_over_k) * lb;

	return height;
}//END cutHeight()

