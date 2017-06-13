/*
 * policy.c
 *
 *  Created on: Apr 13, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "sdlp.h"

extern configType config;

int selectIncumb(incumbType *incumb, omegaType *omega) {

	switch (config.POLICY) {
	case 0:
		return 0;
		break;
	case 1:
		if ( omega == NULL )
			return 0;
		else {
			if ( omega->newPath ) {
				incumb->vals[incumb->cnt] = duplicVector(incumb->vals[0], incumb->len);
				incumb->idx = incumb->cnt++;
				incumb->chg = TRUE;
				return incumb->idx;
			}
			else
				return omega->pathCurrent+1;
		}
		break;
	case 2:
		printf("Not ready yet.\n");
		break;
	default:
		errMsg("algorithm", "selectIncum", "unknown policy type in configuration file", 0);
		break;
	}

	return -1;
}//END selectIncumb()

void checkImprovement(probType **prob, cellType **cell, int numStages) {
	double	candidEst;
	int 	t;

	for ( t = 0; t < numStages-1; t++ ) {
		/* Calculate height at new candidate x with newest cut included */
		candidEst = maxCutHeight(cell[t]->cuts, cell[t]->lb, cell[t]->k, prob[t+1]->coord->colsC, prob[t+1]->num->cntCcols, cell[t]->candidU);
		candidEst += vXvSparse(cell[t]->candidU, prob[t]->dBar);

		/* Calculate height at current incumbent x with newest cut */
		cell[t]->incumb->est[cell[t]->incumb->idx] = maxCutHeight(cell[t]->cuts, cell[t]->lb, cell[t]->k, prob[t+1]->coord->colsC,
				prob[t+1]->num->cntCcols, cell[t]->incumb->vals[0]);
		cell[t]->incumb->est[cell[t]->incumb->idx] += vXvSparse(cell[t]->incumb->vals[0], prob[t]->dBar);

#if VERBOSE
		if ( t == 0 )
			printf("\nEstimates Candidate = %lf\tIncumbent = %lf\n", candidEst, cell[t]->incumb->est[0]);
#endif

		if ((candidEst - cell[t]->incumb->est[cell[t]->incumb->idx]) < (config.R1 * cell[t]->incumb->improv)) {
			/* incumbent update is recommended */
			copyVector(cell[t]->candidU, cell[t]->incumb->vals[cell[t]->incumb->idx], prob[t]->num->cols, TRUE);
			cell[t]->incumb->est[cell[t]->incumb->idx] = candidEst;

			if (cell[t]->k > 1 && cell[t]->incumb->normd_k > 1000*config.TOLERANCE)
				if (cell[t]->incumb->normd_k >= config.R3 * cell[t]->incumb->normd_k_1) {
					cell[t]->incumb->quadScalar *= config.R2 * config.R3 * cell[t]->incumb->normd_k_1 / cell[t]->incumb->normd_k;

					cell[t]->incumb->quadScalar = min(config.MAX_QUAD_SCALAR, cell[t]->incumb->quadScalar );
					cell[t]->incumb->quadScalar = max(config.MIN_QUAD_SCALAR, cell[t]->incumb->quadScalar);
				}

			cell[t]->incumb->normd_k_1 = cell[t]->incumb->normd_k;
			cell[t]->incumb->chg = TRUE;
			printf("+"); fflush(stdout);
		}
		else {
			/* update quad_scalar when incumbent is not updated */
			cell[t]->incumb->quadScalar = min(config.MAX_QUAD_SCALAR, cell[t]->incumb->quadScalar / config.R2);
			cell[t]->incumb->normd_k_1 = cell[t]->incumb->normd_k;
			cell[t]->incumb->chg = FALSE;
		}
	}

#ifdef ALGO_CHECK
		printf("Incumbent estimates = ");
		for ( t = 0; t < T-1; t++ )
			printf("%lf\t", cell[t]->incumbEst[cell[t]->incumbIdx]);
		printf("\n");
#endif

}//END checkImprove()

incumbType *newIncumb(int maxIncumb, vector meanSol, int lenX, sparseVector *dBar) {
	incumbType *incumb;

	if ( !(incumb = (incumbType *) mem_malloc(sizeof(incumbType))) )
		errMsg("allocation", "newIncumb", "incumb", 0);
	if ( !(incumb->vals = (vector *) arr_alloc(maxIncumb, vector)) )
		errMsg("allocation", "newIncumb", "incumb->vals", 0);
	if ( !(incumb->est = (vector) arr_alloc(maxIncumb, double)) )
		errMsg("allocation", "newIncumb", "incumb->est", 0);
	if ( !(incumb->cutidx = (intvec) arr_alloc(maxIncumb, int)) )
		errMsg("allocation", "newIncumb", "incumb->cutidx", 0);
	if ( !(incumb->updtIter = (intvec) arr_alloc(maxIncumb, int)) )
		errMsg("allocation", "newIncumb", "incumb->updtIter", 0);
	incumb->quadScalar = config.MIN_QUAD_SCALAR;
	incumb->len = lenX;

	/* initialize mean value solution as the first incumbent */
	incumb->vals[0] = duplicVector(meanSol, lenX);
	incumb->est[0] = vXvSparse(incumb->vals[0], dBar);
	incumb->cnt = 1;
	incumb->chg = TRUE;
	incumb->normd_k = 0.0;
	incumb->normd_k_1 = 0.0;
	incumb->improv = 0.0;
	incumb->idx = 0;
	incumb->cutidx[0] = -1; //TODO

	return incumb;
}//END newIncumb()

void freeIncumb(incumbType *incumb) {
	int n;

	if ( incumb->vals ) {
		for ( n = 0; n < incumb->cnt; n++ )
			if ( incumb->vals[n] ) mem_free(incumb->vals[n]);
		mem_free(incumb->vals);
	}
	if ( incumb->est ) mem_free(incumb->est);
	if ( incumb->cutidx) mem_free(incumb->cutidx);
	if ( incumb->updtIter ) mem_free(incumb->updtIter);
	mem_free(incumb);

}//END freeIncumb()
