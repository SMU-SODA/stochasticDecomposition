/*
 * stocUpdate.c
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

/* This function calculates a new column in the delta structure, based on a new observation of omega. Thus, lambda_pi X C and lambda_pi X b
 * are calculated for all values of lambda_pi, for the new C(omega) and b(omega).  Room in the array has already been allocated, so the function
 * only fills it, in the column specified by _obs_. It is assumed that this observation is distinct from all previous ones, and thus a new column
 * must be calculated. */
void calcDeltaCol(numType *num, coordType *coord, lambdaType *lambda, vector observ, int omegaIdx, deltaType *delta) {
    int piIdx;
    sparseVector bomega;
    sparseMatrix Comega;
    vector lambPi;
    vector piCrossC;

    bomega.cnt = num->rvbOmCnt;	bomega.col = coord->omegaRow; bomega.val= observ;
    Comega.cnt = num->rvCOmCnt; Comega.col = coord->omegaCol + num->rvbOmCnt;
    Comega.row = coord->omegaRow + num->rvbOmCnt; Comega.val = observ + num->rvbOmCnt;

    /* For all dual vectors, lambda(pi), calculate pi X bomega and pi X Comega */
    for (piIdx = 0; piIdx < lambda->cnt; piIdx++) {
        /* Retrieve a new (sparse) dual vector, and expand it into a full vector */
        lambPi = expandVector(lambda->vals[piIdx], coord->rvRows, num->rvRowCnt, num->rows);

        /* Multiply the dual vector by the observation of bomega and Comega */
        /* Reduce PIxb from its full vector form into a sparse vector */
        delta->vals[piIdx][omegaIdx].pib = vXvSparse(lambPi, &bomega);
        if ( num->rvColCnt != 0 ) {
        	piCrossC = vxMSparse(lambPi, &Comega, num->prevCols);
        	delta->vals[piIdx][omegaIdx].piC = reduceVector(piCrossC, coord->rvCols, num->rvColCnt);
            mem_free(piCrossC);
        }
        else
        	delta->vals[piIdx][omegaIdx].piC = NULL;

        mem_free(lambPi);
    }

}//END calcDeltaCol

/* This function stores a new lambda_pi vector in the lambda structure.  Each lambda_pi represents only those dual variables whose rows in the
 * constraint matrix have random elements.  Thus  the (full) dual vector, Pi,  passed to the function is converted into the sparse vector lambda_pi.
 * This vector is then compared with all previous lambda_pi vectors, searching for a duplication. If a duplicate is found, the vector is not added
 * to the structure, and the function returns the index of the duplicate vector. Otherwise, it adds the vector to the end of the structure,
 *and returns an index to the last element in lambda. */
int calcLambda(numType *num, coordType *coord, vector Pi, lambdaType *lambda, BOOL *newLambdaFlag) {
    int 	pi_idx;
    vector	lambda_pi;

    /* Pull out only those elements in dual vector which have rv's */
    lambda_pi = reduceVector(Pi, coord->rvRows, num->rvRowCnt);

    /* Compare resulting lambda_pi with all previous vectors */
    for (pi_idx = 0; pi_idx < lambda->cnt; pi_idx++)
        if (equalVector(lambda_pi, lambda->vals[pi_idx], num->rvRowCnt, config.TOLERANCE)) {
            mem_free(lambda_pi);
            *newLambdaFlag = FALSE;
            return pi_idx;
        }

    /* Add the vector to lambda structure */
    lambda->vals[lambda->cnt] = lambda_pi;
    *newLambdaFlag = TRUE;

    return lambda->cnt++;
}//END calcLambda

int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, vector pi, double mubBar,
              int idxLambda, BOOL newLambdaFlag, int iter, sigmaType *sigma, BOOL *newSigmaFlag) {
    vector	piCBar, temp;
    double 	pibBar;
    int 	cnt;

    /* sigma = \pi_t^\top \bar{b}_t - \bar{C}_t^\top \pi_t */
    pibBar = vXvSparse(pi, bBar) + mubBar;

    temp = vxMSparse(pi, CBar, num->prevCols);
    piCBar = reduceVector(temp, coord->colsC, num->cntCcols);
    mem_free(temp);

    if (!newLambdaFlag){
        for (cnt = 0; cnt < sigma->cnt; cnt++) {
            if (DBL_ABS(pibBar - sigma->vals[cnt].pib) <= config.TOLERANCE) {
                if (equalVector(piCBar, sigma->vals[cnt].piC, num->cntCcols, config.TOLERANCE))
                    if(sigma->lambdaIdx[cnt]== idxLambda){
                        mem_free(piCBar);
                        (*newSigmaFlag) = FALSE;
                        return cnt;
                    }
            }
        }
    }

    (*newSigmaFlag) = TRUE;
    sigma->vals[sigma->cnt].pib  = pibBar;
    sigma->vals[sigma->cnt].piC  = piCBar;
    sigma->lambdaIdx[sigma->cnt] = idxLambda;
    sigma->ck[sigma->cnt] = iter;

    return sigma->cnt++;

}//END calcSigma()

/* This function calculates a new row in the delta structure, based on a new dual vector, lambda_pi, by calculating lambda_pi X b and
 * lambda_pi X C for all previous realizations of b(omega) and C(omega).  It is assumed that the lambda vector is distinct from all previous ones
 * and thus a new row is warranted. */
int calcDeltaRow(int maxIter, numType *num, coordType *coord, omegaType *omega, lambdaType *lambda, int lambdaIdx, deltaType *delta) {
    sparseVector bomega;
    sparseMatrix Comega;
    vector 	lamb_pi, pixC;
    int		obs;

    bomega.cnt = num->rvbOmCnt;	bomega.col = coord->omegaRow;
    Comega.cnt = num->rvCOmCnt; Comega.col = coord->omegaCol + num->rvbOmCnt; Comega.row = coord->omegaRow + num->rvbOmCnt;

    if ( !(delta->vals[lambdaIdx] = (pixbCType *) arr_alloc(maxIter, pixbCType)))
        errMsg("allocation", "calcDeltaRow", "delta->val[cnt]", 0);

    /* expand the compressed lambda vector */
    lamb_pi = expandVector(lambda->vals[lambdaIdx], coord->rvRows, num->rvRowCnt, num->rows);

    /* go through all the observations and compute pi x b and pi x C */
    for (obs = 0; obs < omega->cnt; obs++) {

        bomega.val= omega->vals[obs];
        Comega.val = omega->vals[obs] + num->rvbOmCnt;

        delta->vals[lambdaIdx][obs].pib = vXvSparse(lamb_pi, &bomega);
        if ( num->rvColCnt != 0 ) {
        	pixC = vxMSparse(lamb_pi, &Comega, num->prevCols);
        	delta->vals[lambdaIdx][obs].piC = reduceVector(pixC, coord->rvCols, num->rvColCnt);
            mem_free(pixC);
        }
        else
        	delta->vals[lambdaIdx][obs].piC = NULL;
    }

    mem_free(lamb_pi);

    return 0;

}//END calcDeltaRow()

/* This function obtains a new vector of realizations of the random variables. It compares the new vector with all previous vectors, looking for
 * a duplication.  If it finds a duplicate, it returns the index of that duplicate; otherwise, it adds the vector to the list of distinct realizations
 * and returns the index of that realization. Note that the simulated observation does not have contain one-norm, while the values stored in
 * omegaType do */
int calcOmega(vector observ, int begin, int end, omegaType *omega, BOOL *newOmegaFlag) {
    int cnt;

    /* Compare vector with all the previous observations */
    for (cnt = 0; cnt < omega->cnt; cnt++)
        if (equalVector(observ, omega->vals[cnt], end-begin, config.TOLERANCE)) {
            (*newOmegaFlag) = FALSE;
            omega->weight[cnt]++;
            return cnt;
        }

    /* Add the realization vector to the list */
    omega->vals[omega->cnt] = duplicVector(observ, end-begin);
    omega->weight[omega->cnt] = 1;
    (*newOmegaFlag) = TRUE;

#ifdef STOCH_CHECK
    printf("Observation (%d): ", *newOmegaFlag);
    printVector(omega->vals[omega->cnt], end - begin, NULL);
#endif

    return omega->cnt++;
}//calcOmega()

/* This function compute the reduced cost of every second stage variables. They will be used to calculate the \mu x b and then added to the \pi x b. */
int computeMU(LPptr lp, intvec cstat, int numCols, double *mubBar) {
    vector	dj, u;
    int		n;

    (*mubBar) = 0.0;

    if ( !(dj = (vector) arr_alloc(numCols+1, double)))
        errMsg("allocation", "computeMu", "dual slacks", 0);
    if ( !(u = (vector) arr_alloc(numCols+1, double)))
        errMsg("allocation", "computeMu", "TDA solutions", 0);

    if ( getPrimal(lp, u, numCols) ) {
        errMsg("solver", "forOptPass", "failed to obtain primal solution", 0);
        return 1;
    }
    if (getDualSlacks(lp, dj, numCols) ) {
        errMsg("solver", "computeMu", "failed to obtain dual slacks", 0);
        return 1;
    }

    for (n = 1; n <= numCols;  n++) {
        switch (cstat[n]) {
            case AT_LOWER:
                (*mubBar) += dj[n]*u[n];
                break;
            case AT_UPPER:
                (*mubBar) += dj[n]*u[n];
                break;
            default:
                break;
        }
    }

    mem_free(u); mem_free(dj);

    return 0;
}//END compute_mu()

oneBasis *newBasis(int maxPhiLength, unsigned long *codedCol, unsigned long *codedRow, int numCols) {
	oneBasis *B;

	if ( !(B = (oneBasis *) mem_malloc(sizeof(oneBasis))))
		errMsg("allocation", "newBasis", "B", 0);
	if ( maxPhiLength > 0 ) {
		if ( !(B->phiHeader = (intvec) arr_alloc(maxPhiLength+1, int)) )
			errMsg("allocation", "newBasis", "B->phiHeader", 0);
		if ( !(B->phi = (vector *) arr_alloc(maxPhiLength, vector)) )
			errMsg("allocation", "newBasis", "B->phi", 0);
		if ( !(B->omegaIdx = (intvec) arr_alloc(maxPhiLength+1, int)) )
			errMsg("allocation", "newBasis", "B->omegaIdx", 0);
		if ( !(B->g = (vector) arr_alloc(numCols+1, double)) )
			errMsg("allocation", "newBasis", "B->g", 0);
		B->cCode  = codedCol;
		B->rCode  = codedRow;
	}
	else {
		B->phiHeader = NULL; B->phi   = NULL;
		B->cCode 	 = NULL; B->rCode = NULL;
		B->psi 		 = NULL;
	}
	if ( !(B->lambdaIdx = (intvec) arr_alloc(maxPhiLength+1, int)) )
		errMsg("allocation", "newBasis", "B->lambdaIdx", 0);
	if ( !(B->sigmaIdx = (intvec) arr_alloc(maxPhiLength+1, int)) )
		errMsg("allocation", "newBasis", "B->sigmaIdx", 0);
	B->weight = 1;
	B->phiLength = 0;

	return B;
}//END newBasis()

/* This function allocates a new basisType data structure which holds all the unique basis discovered by the algorithm. It returns a pointer to the
 * structure. */
basisType *newBasisType(int numIter, int numCols, int numRows, int wordLength) {
	basisType *basis;

	if ( !(basis = (basisType *) mem_malloc(sizeof(basisType))))
		errMsg("allocation", "newBasisType", "basis", 0);
	if ( !(basis->vals = (oneBasis **) arr_alloc(numIter, oneBasis *)))
		errMsg("allocation", "newBasisType", "basis->vals", 0);
	if ( !(basis->obsFeasible = (BOOL **) arr_alloc(numIter, BOOL *)))
		errMsg("allocation", "newBasisType", "basis->obsFeasible", 0);
	basis->cnt = 0;
	basis->cCodeLen = ceil(numCols/wordLength) + 1;
	basis->rCodeLen = ceil(numRows/wordLength) + 1;

	return basis;
}//END newBasis()

void freeBasisType(basisType *basis) {
	int n;

	if ( basis ) {
		if ( basis->vals ) {
			for ( n = 0; n < basis->cnt; n++ ) {
				freeOneBasis(basis->vals[n]);
				mem_free(basis->obsFeasible[n]);
			}
			mem_free(basis->vals);
			mem_free(basis->obsFeasible);
		}
		mem_free(basis);
	}

}//END freeBasisType

void freeOneBasis(oneBasis *B) {
	int n;

	if ( B ) {
		if (B->cCode) mem_free(B->cCode);
		if (B->rCode) mem_free(B->rCode);
		if (B->lambdaIdx) mem_free(B->lambdaIdx);
		if (B->sigmaIdx) mem_free(B->sigmaIdx);
		if (B->omegaIdx) mem_free(B->omegaIdx);
		if (B->phiHeader) mem_free(B->phiHeader);
		if (B->g) mem_free(B->g);
		if ( B->phi) {
			for ( n = 0; n < B->phiLength; n++ )
				if ( B->phi[n] ) mem_free(B->phi[n]);
			mem_free(B->phi);
		}
		if ( B->psi ) freeSparseMatrix(B->psi);
		mem_free(B);
	}

}//END freeOneBasis

/* This function allocates a new lambda structure, with room for num_lambdas lambda vectors of size vect_size.  It returns a pointer to the structure.
 * Only some of the individual lambda vectors are expected to be allocated (according to the num_vect parameter) so that there is room for new
 * lambdas to be created. */
lambdaType *newLambda(int maxLambda, int numLambda, int numRVrows) {
    lambdaType *lambda;
    int cnt;

    if (!(lambda = (lambdaType *) mem_malloc (sizeof(lambdaType))))
        errMsg("allocation", "new_lambda", "lambda",0);

    if (!(lambda->vals = arr_alloc(maxLambda, vector)))
        errMsg("allocation", "new_lambda", "lambda->val",0);

    for (cnt = 0; cnt < numLambda; cnt++)
        if (!(lambda->vals[cnt] = arr_alloc(numRVrows + 1, double)))
            errMsg("allocation", "new_lambda", "lambda->val[cnt]",0);

    lambda->cnt = numLambda;

    return lambda;
}//END new_lambda

/* This function creates a new sigma structure, and allocates memory for the arrays associated with it.  It returns a pointer to this structure.
 * Some pi X T vectors are also allocated, according to the num_vals parameter  (num_vals is expected to be less than num_sigmas, so that there
 * is room for further work).  Note that  memory for sigma->col is not allocated, but is taken from prob.*/
sigmaType *newSigma(int numIter, int numNzCols, int numPi) {
    sigmaType *sigma;
    int cnt;

    if (!(sigma = (sigmaType *) mem_malloc (sizeof(sigmaType))))
        errMsg("allocation", "new_sigma", "sigma",0);
    if (!(sigma->lambdaIdx = (intvec) arr_alloc(numIter, int)))
        errMsg("allocation", "new_sigma", "sigma->lambIdx",0);
    if (!(sigma->ck = (intvec) arr_alloc(numIter, int)))
        errMsg("allocation", "new_sigma", "sigma->ck",0);
    if (!(sigma->vals = arr_alloc(numIter, pixbCType)))
        errMsg("allocation", "new_sigma", "sigma->vals",0);
    for (cnt = 0; cnt < numPi && cnt < numIter; cnt++)
        if (!(sigma->vals[cnt].piC = arr_alloc(numNzCols+1, double)))
            errMsg("allocation", "new_sigma", "sigma->val[cnt]",0);

    sigma->cnt = numPi;

    return sigma;
}//END newSigma

/***********************************************************************\
 ** This function creates a new delta structure with arrays of the specified
 ** size and returns a pointer to it.  Note that the pi X T vectors
 ** themselves are not allocated, since they will not all be filled with
 ** values.  (they are only filled as they are produced).
 ** Not even the arrays of pi_R_T_types are allocated, as this also
 ** occurs in calc_delta_row().  However, the column coordinates of the
 ** eventual multiplications are initialized, since they are known.
 \***********************************************************************/
deltaType *newDelta(int numIter) {
    deltaType *delta;

    if (!(delta = (deltaType *) mem_malloc (sizeof(deltaType))))
        errMsg("Allocation", "new_delta", "d",0);
    if (!(delta->vals = (pixbCType **) arr_alloc(numIter, pixbCType *)))
        errMsg("Allocation", "new_delta", "d->val",0);
    return delta;
}//END newDelta

/* This function allocates memory for an omega structure.  It allocates the memory to structure elements: a vector to hold an array of
 * observation and the weights associated with it. */
omegaType *newOmega(int numIter) {
    omegaType *omega;

    if ( !(omega = (omegaType *) mem_malloc(sizeof(omegaType))) )
        errMsg("allocation","newOmega", "omega", 0);
    if ( !(omega->weight = (intvec) arr_alloc(numIter, int)) )
        errMsg("allocation", "newOmega", "omega->weight", 0);
    if ( !(omega->vals = (vector *) arr_alloc(numIter, vector)) )
        errMsg("allocation", "newOmega", "omega->vals", 0);
    omega->cnt = 0;

    return omega;
}//END newOmega()

void freeOmegaType(omegaType *omega) {
	int n;

	if ( omega->weight ) mem_free(omega->weight);
	if ( omega->vals ) {
		for ( n = 0; n < omega->cnt; n++ )
			if ( omega->vals[n] ) mem_free(omega->vals[n]);
		mem_free(omega->vals);
	}
	mem_free(omega);

}//END freeOmegaType()

void freeLambdaType(lambdaType *lambda) {
	int n;

	if (lambda) {
		if (lambda->vals) {
			for ( n = 0; n < lambda->cnt; n++ )
				if (lambda->vals[n]) mem_free(lambda->vals[n]);
			mem_free(lambda->vals);
		}
		mem_free(lambda);
	}

}//END freeLambdaType()

void freeSigmaType(sigmaType *sigma) {
	int n;

	if (sigma) {
		if (sigma->lambdaIdx) mem_free(sigma->lambdaIdx);
		for ( n = 0; n < sigma->cnt; n++ )
			if (sigma->vals[n].piC) mem_free(sigma->vals[n].piC);
		if (sigma->vals) mem_free(sigma->vals);
		if (sigma->ck) mem_free(sigma->ck);
		mem_free(sigma);
	}

}//END freeSigmaType()

void freeDeltaType (deltaType *delta, int lambdaCnt, int omegaCnt) {
	int n, m;

	if (delta) {
		if (delta->vals) {
			for ( n = 0; n < lambdaCnt; n++ ) {
				if (delta->vals[n]) {
					for ( m = 0; m < omegaCnt; m++ )
						if (delta->vals[n][m].piC)
							mem_free(delta->vals[n][m].piC);
					mem_free(delta->vals[n]);
				}
			}
			mem_free(delta->vals);
		}
		mem_free(delta);
	}

}//END freeDeltaType()
