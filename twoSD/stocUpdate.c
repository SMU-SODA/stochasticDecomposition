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

/* This function updates all the structures necessary for forming a stochastic cut. The latest observation of omega and the latest dual solution
 * to the subproblem are added to their appropriate structures. Then Pi x b and Pi x C are computed and for the latest omega and dual vector,
 * and are added to the appropriate structures.
 * Note that the new column of delta is computed before a new row in lambda is calculated and before the new row in delta is completed,
 * so that the intersection of the new row and new column in delta is only computed once (they overlap at the bottom, right-hand corner). */
int stochasticUpdates(probType *prob, LPptr spLP, lambdaType *lambda, sigmaType *sigma, deltaType *delta, int deltaRowLength, omegaType *omega,
		int omegaIdx, BOOL newOmegaFlag, int currentIter, double TOLERANCE) {
	vector 	piS;
	intvec	cstat;
	double	mubBar;
    int 	lambdaIdx, sigmaIdx;
    BOOL 	newLambdaFlag= FALSE, newSigmaFlag= FALSE;

	/* Allocate memory. */
	if ( !(cstat = (intvec) arr_alloc( prob->num->cols+1, int)))
		errMsg("allocation", "stochasticUpdates", "cstat", 0);
	if ( !(piS = (vector) arr_alloc(prob->num->rows+1, double)) )
		errMsg("allocation", "stochasticUpdates", "piDet", 0);

	/* Obtain the status of columns and rows in the basis. */
	if ( getBasis(spLP, cstat, NULL) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to get the basis column and row status", 0);
		return -1;
	}
	/* Record the dual and reduced cost on bounds. */
	if ( getDual(spLP, piS, prob->num->rows) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to get the dual", 0);
		return 1;
	}

	if ( computeMU(spLP, cstat,  prob->num->cols, &mubBar) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to compute mubBar for subproblem", 0);
		return -1;
	}

    /* Only need to calculate column if new observation of omega found */
    if (newOmegaFlag)
    	calcDelta(prob->num, prob->coord, lambda, omega, delta, deltaRowLength, omegaIdx, newOmegaFlag);

    /* extract the dual solutions corresponding to rows with random elements in them */
    lambdaIdx = calcLambda(prob->num, prob->coord, piS, lambda, &newLambdaFlag, TOLERANCE);

    /* compute Pi x bBar and Pi x Cbar */
    sigmaIdx = calcSigma(prob->num, prob->coord, prob->bBar, prob->Cbar, piS, mubBar, lambdaIdx, newLambdaFlag,
    		currentIter, sigma, &newSigmaFlag, TOLERANCE);

    /* Only need to calculate row if a distinct lambda was found. We could use Pi, instead of lambda(Pi), for this calculation, */
    /* and save the time for expanding/reducing vector even though the lambda is the same, the current Pi might be a
     distinct one due to the variations in sigma*/
    if (newLambdaFlag)
    	calcDelta(prob->num, prob->coord, lambda, omega, delta, deltaRowLength, lambdaIdx, FALSE);

    mem_free(piS); mem_free(cstat);
    return sigmaIdx;
}//END stochasticUpdates

/*This function loops through all the dual vectors found so far and returns the index of the one which satisfies the expression:
 * 				argmax { Pi x (R - T x X) | all Pi }
 * where X, R, and T are given.  It is calculated in this form:
 * 				Pi x bBar + Pi x bomega + (Pi x Cbar) x X + (Pi x Comega) x X.
 * Since the Pi's are stored in two different structures (sigma and delta), the index to the maximizing Pi is actually a structure
 * containing two indices.  (While both indices point to pieces of the dual vectors, sigma and delta may not be in sync with one
 * another due to elimination of non-distinct or redundant vectors. */
int computeIstar(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, vector piCbarX, vector Xvect, int obs,
		int numSamples, BOOL pi_eval, double *argmax, BOOL isNew) {
	double 	arg;
	int 	cnt, maxCnt, sigmaUp, sigmaLow;

	if (pi_eval == TRUE)
		numSamples -= (numSamples / 10 + 1);

	/* Establish the range of iterations over which the istar calculations are conducted. Only bases discovered in this iteration range are used. */
	if ( !isNew ) {
		sigmaUp = numSamples; sigmaLow = -INT_MAX;
	}
	else {
		sigmaUp = INT_MAX; sigmaLow = numSamples;
	}

	*argmax = -DBL_MAX; maxCnt = 0;
	for (cnt = 0; cnt < sigma->cnt; cnt++) {
		if ( sigma->ck[cnt] > sigmaLow && sigma->ck[cnt] <= sigmaUp ) {
			/* Start with (Pi x bBar) + (Pi x bomega) + (Pi x Cbar) x X */
			arg = sigma->vals[cnt].pib + delta->vals[sigma->lambdaIdx[cnt]][obs].pib - piCbarX[cnt];

			/* Subtract (Pi x Comega) x X. Multiply only non-zero VxT values */
			arg -= vXv(delta->vals[sigma->lambdaIdx[cnt]][obs].piC, Xvect, coord->rvCols, num->rvColCnt);

			if (arg > (*argmax)) {
				*argmax = arg;
				maxCnt = cnt;
			}
		}
	}

	if ( (*argmax == -DBL_MAX ) )
		return -1;
	else
		return maxCnt;
}//END computeIstar

/* This function calculates a new column in the delta structure, based on a new observation of omega. Thus, lambda_pi X C and lambda_pi X b
 * are calculated for all values of lambda_pi, for the new C(omega) and b(omega).  Room in the array has already been allocated, so the function
 * only fills it, in the column specified by _obs_. It is assumed that this observation is distinct from all previous ones, and thus a new column
 * must be calculated. */
void calcDelta(numType *num, coordType *coord, lambdaType *lambda, omegaType *omega, deltaType *delta, int deltaRowLength, int elemIdx,
		BOOL newOmegaFlag) {
    int 	idx;
    sparseVector bomega;
    sparseMatrix Comega;
    vector 	lambdaPi, piCrossC;

    bomega.cnt = num->rvbOmCnt;	bomega.col = coord->rvbOmRows;
    Comega.cnt = num->rvCOmCnt; Comega.col = coord->rvCOmCols + num->rvbOmCnt;
    Comega.row = coord->rvCOmRows + num->rvbOmCnt;

    if ( newOmegaFlag ) {
		/* Case I: New observation encountered. */
        bomega.val= omega->vals[elemIdx];
        Comega.val = omega->vals[elemIdx] + num->rvbOmCnt;

        /* For all dual vectors, lambda(pi), calculate pi X bomega and pi X Comega */
        for (idx = 0; idx < lambda->cnt; idx++) {
            /* Retrieve a new (sparse) dual vector, and expand it into a full vector */
            lambdaPi = expandVector(lambda->vals[idx], coord->rvRows, num->rvRowCnt, num->rows);

            /* Multiply the dual vector by the observation of bomega and Comega */
            /* Reduce PIxb from its full vector form into a sparse vector */
            delta->vals[idx][elemIdx].pib = vXvSparse(lambdaPi, &bomega);
            if ( num->rvColCnt != 0 ) {
            	piCrossC = vxMSparse(lambdaPi, &Comega, num->prevCols);
            	delta->vals[idx][elemIdx].piC = reduceVector(piCrossC, coord->rvCols, num->rvColCnt);
                mem_free(piCrossC);
            }
            else
            	delta->vals[idx][elemIdx].piC = NULL;

            mem_free(lambdaPi);
        }
    }
    else {
		/* Case II: New dual vector encountered. */
        if ( !(delta->vals[elemIdx] = (pixbCType *) arr_alloc(deltaRowLength, pixbCType)))
            errMsg("allocation", "calcDeltaRow", "delta->val[cnt]", 0);

        /* expand the compressed lambda vector */
        lambdaPi = expandVector(lambda->vals[elemIdx], coord->rvRows, num->rvRowCnt, num->rows);

        /* go through all the observations and compute pi x b and pi x C */
        for (idx = 0; idx < omega->cnt; idx++) {

            bomega.val= omega->vals[idx];
            Comega.val = omega->vals[idx] + num->rvbOmCnt;

            delta->vals[elemIdx][idx].pib = vXvSparse(lambdaPi, &bomega);
            if ( num->rvColCnt != 0 ) {
            	piCrossC = vxMSparse(lambdaPi, &Comega, num->prevCols);
            	delta->vals[elemIdx][idx].piC = reduceVector(piCrossC, coord->rvCols, num->rvColCnt);
                mem_free(piCrossC);
            }
            else
            	delta->vals[elemIdx][idx].piC = NULL;
        }
        mem_free(lambdaPi);
    }

}//END calcDelta()

/* This function stores a new lambda_pi vector in the lambda structure.  Each lambda_pi represents only those dual variables whose rows in the
 * constraint matrix have random elements.  Thus  the (full) dual vector, Pi,  passed to the function is converted into the sparse vector lambda_pi.
 * This vector is then compared with all previous lambda_pi vectors, searching for a duplication. If a duplicate is found, the vector is not added
 * to the structure, and the function returns the index of the duplicate vector. Otherwise, it adds the vector to the end of the structure,
 *and returns an index to the last element in lambda. */
int calcLambda(numType *num, coordType *coord, vector Pi, lambdaType *lambda, BOOL *newLambdaFlag, double TOLERANCE) {
    int 	pi_idx;
    vector	lambda_pi;

    /* Pull out only those elements in dual vector which have rv's */
    lambda_pi = reduceVector(Pi, coord->rvRows, num->rvRowCnt);

    /* Compare resulting lambda_pi with all previous vectors */
    for (pi_idx = 0; pi_idx < lambda->cnt; pi_idx++)
        if (equalVector(lambda_pi, lambda->vals[pi_idx], num->rvRowCnt, TOLERANCE)) {
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
              int idxLambda, BOOL newLambdaFlag, int currentIter, sigmaType *sigma, BOOL *newSigmaFlag, double TOLERANCE) {
    vector	piCBar, temp;
    double 	pibBar;
    int 	cnt;

    /* sigma = \pi_t^\top \bar{b}_t - \bar{C}_t^\top \pi_t */
    pibBar = vXvSparse(pi, bBar) + mubBar;

    temp = vxMSparse(pi, CBar, num->prevCols);
    piCBar = reduceVector(temp, coord->CCols, num->cntCcols);
    mem_free(temp);

    if (!newLambdaFlag){
        for (cnt = 0; cnt < sigma->cnt; cnt++) {
            if (DBL_ABS(pibBar - sigma->vals[cnt].pib) <= TOLERANCE) {
                if (equalVector(piCBar, sigma->vals[cnt].piC, num->cntCcols, TOLERANCE))
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
    sigma->ck[sigma->cnt] = currentIter;

    return sigma->cnt++;

}//END calcSigma()

/* This function obtains a new vector of realizations of the random variables. It compares the new vector with all previous vectors, looking for
 * a duplication.  If it finds a duplicate, it returns the index of that duplicate; otherwise, it adds the vector to the list of distinct realizations
 * and returns the index of that realization. Note that the simulated observation does not have contain one-norm, while the values stored in
 * omegaType do */
int calcOmega(vector observ, int begin, int end, omegaType *omega, BOOL *newOmegaFlag, double TOLERANCE) {
    int cnt;

    /* Compare vector with all the previous observations */
    for (cnt = 0; cnt < omega->cnt; cnt++)
        if (equalVector(observ, omega->vals[cnt], end-begin, TOLERANCE)) {
            (*newOmegaFlag) = FALSE;
            omega->weights[cnt]++;
            return cnt;
        }

    /* Add the realization vector to the list */
    omega->vals[omega->cnt] = duplicVector(observ, end-begin);
    omega->weights[omega->cnt] = 1;
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

/* This function allocates a new lambda structure, with room for num_lambdas lambda vectors of size vect_size.  It returns a pointer to the structure.
 * Only some of the individual lambda vectors are expected to be allocated (according to the num_vect parameter) so that there is room for new
 * lambdas to be created. */
lambdaType *newLambda(int num_iter, int numLambda, int numRVrows) {
    lambdaType *lambda;
    int cnt;

    if (!(lambda = (lambdaType *) mem_malloc (sizeof(lambdaType))))
        errMsg("allocation", "newLambda", "lambda",0);

    if (!(lambda->vals = arr_alloc(num_iter, vector)))
        errMsg("allocation", "newLambda", "lambda->val",0);

    for (cnt = 0; cnt < numLambda; cnt++)
        if (!(lambda->vals[cnt] = arr_alloc(numRVrows + 1, double)))
            errMsg("allocation", "newLambda", "lambda->val[cnt]",0);

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
        errMsg("allocation", "newSigma", "sigma",0);
    if (!(sigma->lambdaIdx = (intvec) arr_alloc(numIter, int)))
        errMsg("allocation", "newSigma", "sigma->lambIdx",0);
    if (!(sigma->ck = (intvec) arr_alloc(numIter, int)))
        errMsg("allocation", "newSigma", "sigma->ck",0);
    if (!(sigma->vals = arr_alloc(numIter, pixbCType)))
        errMsg("allocation", "newSigma", "sigma->vals",0);
    for (cnt = 0; cnt < numPi && cnt < numIter; cnt++)
        if (!(sigma->vals[cnt].piC = arr_alloc(numNzCols+1, double)))
            errMsg("allocation", "newSigma", "sigma->val[cnt]",0);

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
        errMsg("Allocation", "newDelta", "d",0);
    if (!(delta->vals = (pixbCType **) arr_alloc(numIter, pixbCType *)))
        errMsg("Allocation", "newDelta", "d->val",0);
    return delta;
}//END newDelta

/* This function allocates memory for an omega structure.  It allocates the memory to structure elements: a vector to hold an array of
 * observation and the weights associated with it. */
omegaType *newOmega(int numIter) {
    omegaType *omega;

    if ( !(omega = (omegaType *) mem_malloc(sizeof(omegaType))) )
        errMsg("allocation","newOmega", "omega", 0);
    if ( !(omega->weights = (intvec) arr_alloc(numIter, int)) )
        errMsg("allocation", "newOmega", "omega->weights", 0);
    if ( !(omega->vals = (vector *) arr_alloc(numIter, vector)) )
        errMsg("allocation", "newOmega", "omega->vals", 0);
    omega->cnt = 0;

    return omega;
}//END newOmega()

void freeOmegaType(omegaType *omega, BOOL partial) {
	int n;

	if ( omega->vals ) {
		for ( n = 0; n < omega->cnt; n++ )
			if ( omega->vals[n] ) mem_free(omega->vals[n]);
		if ( partial ) {
			omega->cnt = 0;
			return;
		}
		mem_free(omega->vals);
	}
	if ( omega->weights ) mem_free(omega->weights);
//	if ( omega->probs ) mem_free(omega->probs);
	mem_free(omega);

}//END freeOmegaType()

void freeLambdaType(lambdaType *lambda, BOOL partial) {
	int n;

	if (lambda) {
		if (lambda->vals) {
			for ( n = 0; n < lambda->cnt; n++ )
				if (lambda->vals[n]) mem_free(lambda->vals[n]);
			if ( partial ) {
				lambda->cnt = 0;
				return;
			}
			mem_free(lambda->vals);
		}
		mem_free(lambda);
	}

}//END freeLambdaType()

void freeSigmaType(sigmaType *sigma, BOOL partial) {
	int n;

	if (sigma) {
		for ( n = 0; n < sigma->cnt; n++ )
			if (sigma->vals[n].piC) mem_free(sigma->vals[n].piC);
		if ( partial ) {
			sigma->cnt = 0;
			return;
		}
		if (sigma->lambdaIdx) mem_free(sigma->lambdaIdx);
		if (sigma->vals) mem_free(sigma->vals);
		if (sigma->ck) mem_free(sigma->ck);
		mem_free(sigma);
	}

}//END freeSigmaType()

void freeDeltaType (deltaType *delta, int numDeltaRows, int omegaCnt, BOOL partial) {
	int n, m;

	if (delta) {
		if (delta->vals) {
			for ( n = 0; n < numDeltaRows; n++ ) {
				if (delta->vals[n]) {
					for ( m = 0; m < omegaCnt; m++ )
						if (delta->vals[n][m].piC)
							mem_free(delta->vals[n][m].piC);
					mem_free(delta->vals[n]);
				}
			}
			if ( partial )
				return;
			mem_free(delta->vals);
		}
		mem_free(delta);
	}

}//END freeDeltaType()
