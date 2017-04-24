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

int selectIncumb(incumbType *incumb) {

	switch (config.POLICY) {
	case 0:
		return 0;
		break;
	case 1:
		printf("Not ready yet.\n");
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
		cell[t]->incumb->est[0] = maxCutHeight(cell[t]->cuts, cell[t]->lb, cell[t]->k, prob[t+1]->coord->colsC, prob[t+1]->num->cntCcols,
				cell[t]->incumb->vals[0]);
		cell[t]->incumb->est[0] += vXvSparse(cell[t]->incumb->vals[0], prob[t]->dBar);

#if 0
		if ( t == 0 )
			printf("Estimates Candidate = %lf\tIncumbent = %lf\n", candidEst, cell[t]->incumb->est[0]);
#endif

		if ((candidEst - cell[t]->incumb->est[0]) < (config.R1 * cell[t]->incumb->improv)) {
			/* incumbent update is recommended */
			copyVector(cell[t]->candidU, cell[t]->incumb->vals[0], prob[t]->num->cols, TRUE);
			cell[t]->incumb->est[0] = candidEst;

			if (cell[t]->k > 1 && cell[t]->incumb->normd_k > config.TOLERANCE)
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
