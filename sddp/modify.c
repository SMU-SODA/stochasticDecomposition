/*
 * modify.c
 *
 *  Created on: Dec 5, 2015
 *      Author: Harsha Gangammanavar
 */

#include <utils.h>
#include <smps.h>
#include <prob.h>
#include <solver.h>
#include <sddp.h>

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
}//END computeRHS()
