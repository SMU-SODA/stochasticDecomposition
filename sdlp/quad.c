/*
 * quad.c
 *
 *  Created on: Apr 12, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "sdlp.h"

int constructQP(LPptr lp, int numCols, double regSigma) {
	vector	qpsepvec;
	int		n, status;
#ifdef TRACE
	trPrint("constructQP", 0);
#endif

	/* allocate memory for regularizing terms */
	if ( !(qpsepvec = (vector) arr_alloc(numCols+1, double)) )
		errMsg("allocation", "constructQP", "regularization term vector", 0);

	/* (i) Add regularization term */
	/* matrix with diagonal elements set to $\sigma/2$, the last term corresponds to $\eta$ */
	for ( n = 0; n < numCols; n++ )
		qpsepvec[n] = regSigma;
	qpsepvec[n] = 0;

	/* Now copy the Q matrix for in the solver. */
	status = copyQPseparable(lp, qpsepvec);
	if ( status ) {
		errMsg("solver", "constructQP", "copy Q matrix into the solver", 0);
		return 1;
	}

	mem_free(qpsepvec);
	return 0;
}//END constructQP()

int changeQPrhs(LPptr lp, intvec betaCols, int betaLen, int numRows, sparseMatrix *Dbar, sparseVector *bBar, cutsType *cuts, vector X, vector rhs,
		int numObs, double lb) {
	vector qpRHS;
	intvec indices;
	int		n;

	if ( !(qpRHS = (vector) arr_alloc(numRows, double)) )
		errMsg("allocation", "computeQPrhs", "qpRHS", 0);
	if ( !(indices = (intvec) arr_alloc(numRows, int)))
		errMsg("allocation", "computeQPrhs", "indices", 0);

	/* change the right-hand side to \bar{b}_t - D_t \hat{u}_t */
	for ( n = 0; n < numRows; n++) {
		indices[n]= n;
		qpRHS[n+1] = rhs[n+1];
	}

	/* change the right-hand side using Dbar */
	qpRHS = MSparsexvSub(Dbar, X, qpRHS);

	/* changing cut right-hand side of stage subproblem */
	for ( n = 0; n < cuts->cnt; n++ )
		qpRHS[cuts->vals[n]->rowNum+1] = cuts->vals[n]->alpha - vXv(cuts->vals[n]->beta, X, betaCols, betaLen) +
		((double) numObs / (double) cuts->vals[n]->numObs - 1) * lb;;

	mem_free(indices);
	mem_free(qpRHS);

	return 0;
}//END changeQPrhs()

int changeQPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector X) {
	intvec 	Cindices;
	vector	lbounds, ubounds;
	string	llu, ulu;
	int		n;

	if (!(Cindices = (intvec) arr_alloc(numCols, int)))
		errMsg ("allocation", "change_bounds", "Cindices", 0);
	if (!(lbounds = (vector) arr_alloc(numCols, double)))
		errMsg ("allocation", "change_bounds", "lbounds", 0);
	if (!(llu= (string) arr_alloc(numCols, char)))
		errMsg ("allocation", "change_bounds", "llu", 0);

	if (!(ubounds = (vector) arr_alloc(numCols, double)))
		errMsg ("allocation", "change_bounds", "ubounds", 0);
	if (!(ulu= (string) arr_alloc(numCols, char)))
		errMsg ("allocation", "change_bounds", "ulu", 0);

	for (n = 0; n < numCols; n++) {
		Cindices[n] = n;
		ubounds[n] = bdu[n]- X[n+1];
		ulu[n] = 'U';
		lbounds[n] = bdl[n] - X[n+1];
		llu[n] = 'L';
	}

	if ( changeBDS (lp, numCols, Cindices, ulu, ubounds) ) {
		errMsg("solver", "updtQPbds", "failed to change the upper bounds with incumbent information", 0);
		return 1;
	}

	if ( changeBDS (lp, numCols, Cindices, llu, lbounds) ) {
		errMsg("solver", "updtQPbds", "failed to change the lower bounds with incumbent information", 0);
		return 1;
	}

	mem_free(Cindices);
	mem_free(lbounds); mem_free(llu);
	mem_free(ubounds); mem_free(ulu);

	return 0;
}//END changeQPbds()
