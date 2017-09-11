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

#include "stoc.h"
#include "twoSD.h"

extern configType config;

int stochasticUpdates(cellType *cell, probType *prob, int omegaIdx, BOOL newOmegaFlag) {
	intvec 	cstat, rstat;
	int 	basisIdx, lambdaIdx, sigmaIdx, cnt, offset;
	BOOL	newBasisFlag, newLambdaFlag, newSigmaFlag;

	/* Allocate memory. */
	if ( !(cstat = (intvec) arr_alloc( prob->num->cols+1, int)))
		errMsg("allocation", "getIndexNumber", "cstat", 0);
	if ( !(rstat = (intvec) arr_alloc( prob->num->rows+1, int)))
		errMsg("allocation", "getIndexNumber", "rstat", 0);

	/* Obtain the status of columns and rows in the basis. */
	if ( getBasis(cell->subprob->lp, cstat+1, rstat+1) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to get the basis column and row status", 0);
		return -1;
	}
	/* Record the dual and reduced cost on bounds. */
	if ( getDual(cell->subprob->lp, cell->piS, prob->num->rows) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to get the dual", 0);
		return -1;
	}
	if ( computeMU(cell->subprob->lp, cstat,  prob->num->cols, &cell->mubBar) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to compute mubBar for subproblem", 0);
		return -1;
	}

	/* Update the column of delta structure if a new observation was encountered, and check the feasibility of existing bases with respect to new observation. */
	if ( newOmegaFlag )
		calcDelta(prob->num, prob->coord, cell->basis, cell->lambda, cell->delta, cell->omega, newOmegaFlag, omegaIdx, config.MAX_ITER);

	if ( prob->num->rvdOmCnt > 0 ) {
		/* The random variables corresponding to cost coefficients are listed at the end of vector, the offset is used to index them */
		offset = prob->num->rvbOmCnt + prob->num->rvCOmCnt;

		/* If the cost-coefficients are random, update the basis structure. */
		basisIdx = calcBasis(cell->subprob->lp, cell->basis, prob->dBar, cstat, prob->num->cols, rstat, prob->num->rows,
				prob->coord->rvCols, prob->num->rvdOmCnt, &newBasisFlag, cell->k);

		if ( cell->basis->vals[basisIdx]->phiLength > 0) {
			/* Decompose the dual solution into deterministic and stochastic components. */
			decomposeDualSolution(cell->basis->vals[basisIdx]->phi, cell->omega->vals[omegaIdx]+offset, cell->basis->vals[basisIdx]->omegaIdx,
					cell->basis->vals[basisIdx]->phiLength, cell->piS, prob->num->rows);
		}

		if ( newBasisFlag ) {
			/* Calculations with respect to deterministic component of the dual solution */
			/* Extract the deterministic component of dual solutions corresponding to rows with random elements in them */
			lambdaIdx = cell->basis->vals[basisIdx]->lambdaIdx[0] = calcLambda(prob->num, prob->coord, cell->piS, cell->lambda, &newLambdaFlag);
			if ( newLambdaFlag )
				if ( !(cell->delta->vals[lambdaIdx] = (pixbCType *) arr_alloc(config.MAX_ITER, pixbCType)))
					errMsg("allocation", "stochasticUpdates", "delta->val[cnt]", 0);

			/* Compute the product of deterministic component of dual solution with deterministic (mean value) right-hand side and transfer matrix. */
			cell->basis->vals[basisIdx]->sigmaIdx[0] = calcSigma(prob->num, prob->coord, prob->bBar, prob->Cbar, cell->piS, cell->mubBar,
					lambdaIdx, newLambdaFlag, cell->k, cell->sigma, &newSigmaFlag);

			/* Calculations with respect to stochastic component of the dual solution */
			for (cnt = 0; cnt < cell->basis->vals[basisIdx]->phiLength; cnt++ ) {
				/* Extract the deterministic component of dual solutions corresponding to rows with random elements in them */
				lambdaIdx = cell->basis->vals[basisIdx]->lambdaIdx[cnt+1] = calcLambda(prob->num, prob->coord, cell->basis->vals[basisIdx]->phi[cnt], cell->lambda, &newLambdaFlag);
				if ( newLambdaFlag )
					if ( !(cell->delta->vals[lambdaIdx] = (pixbCType *) arr_alloc(config.MAX_ITER, pixbCType)))
						errMsg("allocation", "stochasticUpdates", "delta->val[cnt]", 0);

				/* Compute the product of deterministic component of dual solution with deterministic (mean value) right-hand side and transfer matrix. */
				cell->basis->vals[basisIdx]->sigmaIdx[cnt+1] = calcSigma(prob->num, prob->coord, prob->bBar, prob->Cbar, cell->basis->vals[basisIdx]->phi[cnt], cell->mubBar,
						lambdaIdx, newLambdaFlag, cell->k, cell->sigma, &newSigmaFlag);
			}

			/* Establish the feasibility of the new basis with respect to all the observations encountered thus far and compute the corresponding delta elements. */
			calcDelta(prob->num, prob->coord, cell->basis, cell->lambda, cell->delta, cell->omega, FALSE, basisIdx, config.MAX_ITER);
		}
	}
	else {
		/* extract the dual solutions corresponding to rows with random elements in them */
		lambdaIdx = calcLambda(prob->num, prob->coord, cell->piS, cell->lambda, &newLambdaFlag);
		if ( newLambdaFlag )
			if ( !(cell->delta->vals[lambdaIdx] = (pixbCType *) arr_alloc(config.MAX_ITER, pixbCType)))
				errMsg("allocation", "stochasticUpdates", "delta->val[cnt]", 0);

		/* compute Pi x bBar and Pi x Cbar */
		sigmaIdx = calcSigma(prob->num, prob->coord, prob->bBar, prob->Cbar, cell->piS, cell->mubBar, lambdaIdx, newLambdaFlag, cell->k, cell->sigma, &newSigmaFlag);

		if ( newSigmaFlag ) {
			basisIdx = cell->basis->cnt++;
			cell->basis->vals[basisIdx] = newBasis(NULL, NULL, NULL, NULL, 0, 0, 0, cell->k, NULL);
			cell->basis->vals[basisIdx]->lambdaIdx[0] = lambdaIdx;
			cell->basis->vals[basisIdx]->sigmaIdx[0]  = sigmaIdx;

			calcDelta(prob->num, prob->coord, cell->basis, cell->lambda, cell->delta, cell->omega, FALSE, basisIdx, config.MAX_ITER);
		}
		else
			basisIdx = sigmaIdx;
	}

	mem_free(cstat);
	mem_free(rstat);
	return basisIdx;
}//End stochasticUpdates()

int calcBasis(LPptr lp, basisType *basis, sparseVector *dBar, intvec cstat, int numCols, intvec rstat, int numRows, intvec rvCols, int rvdOmCnt, BOOL *newBasisFlag, int currentIter) {
	unsigned long *codedCol, *codedRow;
	int		cnt, i;

	/* encode the row and column status */
	codedCol = encodeIntvec(cstat, numCols, WORDLENGTH);
	codedRow = encodeIntvec(rstat, numRows, WORDLENGTH);

	/* check to see if the current basis was encountered before */
	for ( cnt = 0; cnt < basis->cnt; cnt++ ) {
		if ( equalLongIntvec(codedCol, basis->vals[cnt]->cCode, basis->cCodeLen) && equalLongIntvec(codedRow, basis->vals[cnt]->rCode, basis->rCodeLen) ) {
			/* The basis is the same as one encountered before */
			basis->vals[cnt]->weight++;
			mem_free(codedRow); mem_free(codedCol);
#if defined (STOCH_CHECK)
			printf("An old basis encountered :: %d\n", cnt);
#endif
			(*newBasisFlag) = FALSE;
			return cnt;
		}
	}

	/* New basis encountered, add it to the list */
	(*newBasisFlag) = TRUE;
	basis->vals[cnt] = newBasis(lp, codedCol, codedRow, rvCols, numCols, numRows, rvdOmCnt, currentIter, dBar);

#if defined (STOCH_CHECK)
	printf("New basis identified     :: %d\n", cnt);
	printf("\tNumber of basic columns with random cost coefficients      = %d\n", basis->vals[cnt]->phiLength);
	if ( basis->vals[cnt]->phiLength > 0 ) {
		printf("\tIndex in observation vector corresponding to basic columns = "); printIntvec(basis->vals[cnt]->omegaIdx, basis->vals[cnt]->phiLength, NULL);
		printf("\tPhi = ");
		for (i = 0; i < basis->vals[cnt]->phiLength; i++ ) {
			printf("\t\t"); printVector(basis->vals[cnt]->phi[i], numRows, NULL);
		}
	}
	else {
		printf("\tBasic stochastic columns                                   = NULL\n");
		printf("\tPhi = NULL\n");
	}
#endif

	return basis->cnt++;
}//END calcBasis()

int calcDelta(numType *num, coordType *coord, basisType *basis, lambdaType *lambda, deltaType *delta, omegaType *omega, BOOL newOmegaFlag, int elemIdx, int maxIter) {
	int 		 cnt, c, offset[3], lambdaIdx;
	vector 		 lambdaPi, piCrossC;
	sparseVector bOmega, dOmega;
	sparseMatrix COmega;

	/* extract the coordinates and number of random elements */
	bOmega.cnt = num->rvbOmCnt;	bOmega.col = coord->omegaRow;
	offset[1] = num->rvbOmCnt;

	COmega.cnt = num->rvCOmCnt; COmega.col = coord->omegaCol + offset[1]; COmega.row = coord->omegaRow + offset[1];
	offset[2] = offset[1] + num->rvCOmCnt;

	dOmega.cnt = num->rvdOmCnt; dOmega.col = coord->omegaCol+offset[2];

	if ( newOmegaFlag ) {
		/* Case I: New observation encountered. */
		bOmega.val = omega->vals[elemIdx];
		COmega.val = omega->vals[elemIdx] + offset[1];
		dOmega.val = omega->vals[elemIdx] + offset[2];

		/* Loop though all the basis to establish feasibility with respect to new observations, and if feasible, compute delta elements for corresponding lambdas. */
		for ( cnt = 0; cnt < basis->cnt; cnt++ ) {
			/* Establish if the basis is feasible or not. */
			basis->obsFeasible[cnt][elemIdx] = checkBasisFeasibility(basis->vals[cnt], basis->feasSenx, dOmega.val, dOmega.col, dOmega.cnt, num->cols);

			if ( basis->obsFeasible[cnt][elemIdx] ) {
				/* If the basis is feasible, then compute the delta elements. */
				for ( c = 0; c <= basis->vals[cnt]->phiLength; c++ ) {
					lambdaIdx = basis->vals[cnt]->lambdaIdx[c];
					/* Retrieve a new (sparse) dual vector, and expand it into a full vector */
					lambdaPi = expandVector(lambda->vals[lambdaIdx], coord->rvRows, num->rvRowCnt, num->rows);

					/* Multiply the dual vector by the observation of bomega and Comega */
					/* Reduce PIxb from its full vector form into a sparse vector */
					delta->vals[lambdaIdx][elemIdx].pib = vXvSparse(lambdaPi, &bOmega);
					if ( num->rvCOmCnt != 0 ) {
						piCrossC = vxMSparse(lambdaPi, &COmega, num->prevCols);
						delta->vals[lambdaIdx][elemIdx].piC = reduceVector(piCrossC, coord->rvCols, num->rvColCnt);
						mem_free(piCrossC);
					}
					else
						delta->vals[lambdaIdx][elemIdx].piC = NULL;

					mem_free(lambdaPi);
				}
			}
			else {
				for ( c = 0; c <= basis->vals[cnt]->phiLength; c++ ) {
					lambdaIdx = basis->vals[cnt]->lambdaIdx[c];
					delta->vals[lambdaIdx][elemIdx].piC = NULL;
				}
			}
		}
	}
	else {
		/* Case II: New basis encountered. */
		if ( !(basis->obsFeasible[elemIdx] = (BOOL*) arr_alloc(maxIter, BOOL)) )
			errMsg("allocation", "calcDelta", "basis->obsFeasibility[n]", 0);

		/* Loop through all the observations and establish feasibility of new basis with respect to each. If the new basis is feasible compute the delta elements
		 for all lambdas associated with the new basis. */
		for ( cnt = 0; cnt < omega->cnt; cnt++ ) {
			dOmega.val = omega->vals[cnt] + offset[2];

			/* Establish the feasibility of new basis with respect to existing observations */
			basis->obsFeasible[elemIdx][cnt] = checkBasisFeasibility(basis->vals[elemIdx], basis->feasSenx, dOmega.val, dOmega.col, dOmega.cnt, num->cols);

			if ( basis->obsFeasible[elemIdx][cnt] ) {
				/* If the basis is feasible, compute the delta elements */
				bOmega.val = omega->vals[cnt];
				COmega.val = omega->vals[cnt] + offset[1];

				for ( c = 0; c <= basis->vals[elemIdx]->phiLength; c++ ) {
					/* Retrieve a new (sparse) dual vector, and expand it into a full vector */
					lambdaIdx = basis->vals[elemIdx]->lambdaIdx[c];
					lambdaPi = expandVector(lambda->vals[lambdaIdx], coord->rvRows, num->rvRowCnt, num->rows);

					delta->vals[lambdaIdx][cnt].pib = vXvSparse(lambdaPi, &bOmega);
					if ( num->rvCOmCnt != 0 ) {
						piCrossC = vxMSparse(lambdaPi, &COmega, num->prevCols);
						delta->vals[lambdaIdx][cnt].piC = reduceVector(piCrossC, coord->rvCols, num->rvColCnt);
						mem_free(piCrossC);
					}
					else
						delta->vals[lambdaIdx][cnt].piC = NULL;
					mem_free(lambdaPi);
				}
			}
			else {
				for ( c = 0; c <= basis->vals[elemIdx]->phiLength; c++ ) {
					lambdaIdx = basis->vals[elemIdx]->lambdaIdx[c];
					delta->vals[lambdaIdx][cnt].piC = NULL;
				}
			}
		}
	}

	return 0;
}//END calcDelta()

BOOL checkBasisFeasibility(oneBasis *B, vector senx, vector dOmega, intvec rvCols, int rvdOmCnt, int numCols) {
	vector 	reducedCost;
	int 	c;

	if ( !(reducedCost = (vector) arr_alloc(numCols+1, double)) )
		errMsg("allocation", "calcDelta", "costVector", 0);

	if ( rvdOmCnt > 0 ) {
		copyVector(B->gBar, reducedCost, numCols, TRUE);
		addVectors(reducedCost, dOmega, rvCols, rvdOmCnt);
		if ( B->phiLength > 0 ) {
			MSparsexvSub(B->psi, dOmega, reducedCost);
		}
		c = 1;
		while ( c <= numCols) {
			if ( (0*reducedCost[c]) < 0 ) {
				mem_free(reducedCost);
				return FALSE;
			}
			c++;
		}
	}

	mem_free(reducedCost);
	return TRUE;
}//END checkBasisFeasibility()

/* This function calculates a new column in the delta structure, based on a new observation of basis. Thus, lambda_pi X C and lambda_pi X b
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
		if ( num->rvCOmCnt != 0 ) {
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
					if(sigma->lambdaIdx[cnt] == idxLambda){
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
		if ( num->rvCOmCnt != 0 ) {
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

int decomposeDualSolution(vector *phi, vector omegaVals, intvec phiOmegaIdx, int phiLength, vector Pi, int numRows) {
	int n, i;

	for ( n = 0; n < phiLength; n++ )
		for ( i = 1; i <= numRows; i++ )
			Pi[i] -= phi[n][i]*omegaVals[phiOmegaIdx[n+1]];

	return 0;
}//END decomposeDualSolution()

oneBasis *newBasis(LPptr lp, unsigned long *codedCol, unsigned long *codedRow, intvec rvCols, int numCols, int numRows, int rvdOmCnt, int currentIter, sparseVector *dBar) {
	oneBasis *B;
	vector 	 tempPsiRow, costVector, basicCost;
	intvec	 basisHead, phiHead;
	int 	 i, j;

	/* allocate memory to elements of the basis structure */
	if ( !(B = (oneBasis *) mem_malloc(sizeof(oneBasis))))
		errMsg("allocation", "newBasis", "B", 0);
	if ( !(B->lambdaIdx = (intvec) arr_alloc(rvdOmCnt+1, int)) )
		errMsg("allocation", "newBasis", "B->lambdaIdx", 0);
	if ( !(B->sigmaIdx = (intvec) arr_alloc(rvdOmCnt+1, int)) )
		errMsg("allocation", "newBasis", "B->sigmaIdx", 0);
	if ( !(B->omegaIdx = (intvec) arr_alloc(rvdOmCnt+1, int)) )
		errMsg("allocation", "newBasis", "B->lambdaIdx", 0);
	B->cCode  	 = codedCol;
	B->rCode     = codedRow;
	B->ck    	 = currentIter;
	B->weight 	 = 1;
	B->phiLength = 0;
	B->phi 		 = NULL; B->psi = NULL; B->gBar = NULL;

	if ( rvdOmCnt == 0 )
		return B;

	/* Allocate memory for the basis header. */
	if ( !(basisHead = (intvec) arr_alloc(numRows+1, int)) )
		errMsg("allocation", "calcBasis", "basisHead", 0);
	if ( !(phiHead = (intvec) arr_alloc(rvdOmCnt+1, int)) )
		errMsg("allocation", "calcBasis", "basisHead", 0);
	if ( !(B->gBar = (vector) arr_alloc(numCols+1, double)) )
		errMsg("allocation", "newBasis", "B->gBar", 0);

	/* Compute the phi matrix associated with the current basis. We begin by first identifying the basis header. A negative value in basis header indicates a slack row. */
	getBasisHead(lp, basisHead+1, NULL);

	/* Compute the phi matrix header and extract the basis (of the dual) inverse matrix rows corresponding to the header. */
	for ( i = 1; i <= rvdOmCnt; i++ )		/* Loop through all the columns with random cost coefficients to see if any of them are basic */
		if ( (phiHead[i] = isElementIntvec(basisHead, numRows, (rvCols[i]-1) )) > 0 ) {
			/* _idx_ > 0 gives the index of column with random cost coefficient in the basis head */
			if ( B->phiLength == 0 ) {
				if ( !(B->phi = (vector *) arr_alloc(rvdOmCnt, vector)) )
					errMsg("allocation", "newBasis", "B->phi", 0);
			}
			if ( !(B->phi[B->phiLength] = (vector) arr_alloc(numRows+1, double)) )
				errMsg("allocation", "calcBasis", "B->phi[i]", 0);
			getBasisInvRow(lp, phiHead[i]-1, B->phi[B->phiLength]+1);
			B->phiLength++;
			B->omegaIdx[B->phiLength] = i;
		}

	if ( B->phiLength > 0 ) {
		/* Reallocate memory to elements already assigned and initialize for the remainder of the elements. */
		B->phi = (vector *) mem_realloc(B->phi, B->phiLength*sizeof(vector));
		B->omegaIdx = (intvec) mem_realloc(B->omegaIdx, (B->phiLength+1)*sizeof(int));
		B->sigmaIdx = (intvec) mem_realloc(B->sigmaIdx, (B->phiLength+1)*sizeof(int));
		B->lambdaIdx = (intvec) mem_realloc(B->lambdaIdx, (B->phiLength+1)*sizeof(int));

		if ( !(B->psi = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))) )
			errMsg("allocation", "newBasis", "B->psi", 0);
		B->psi->val = (vector) arr_alloc(numCols*B->phiLength+1, double);
		B->psi->col = (intvec) arr_alloc(numCols*B->phiLength+1, int);
		B->psi->row = (intvec) arr_alloc(numCols*B->phiLength+1, int);
		B->psi->cnt = 0;
	}

	/* Extract the basic variable cost vector and psi matrix (the tableau entries) */
	if ( !(basicCost = (vector) arr_alloc(numRows+1, double)) )
		errMsg("allocation", "newBasis", "basicCost", 0);
	costVector = expandVector(dBar->val, dBar->col, dBar->cnt, numCols);
	for ( i = 1; i <= numRows; i++ )
		basicCost[i] = costVector[basisHead[i]+1];

	if ( !(tempPsiRow = (vector) arr_alloc(numRows+1, double)) )
		errMsg("allocation", "newBasis", "tempPsiRow", 0);

	for ( i = 1; i <= numCols; i++ ) {
		getBasisInvACol(lp, i-1, tempPsiRow+1);

		B->gBar[i] = costVector[i] - vXv(tempPsiRow, basicCost, NULL, numRows);

		for ( j = 1; j <= B->phiLength; j++ ) {
			B->psi->row[B->psi->cnt+1] = i;
			B->psi->col[B->psi->cnt+1] = B->omegaIdx[j];
			B->psi->val[B->psi->cnt+1] = tempPsiRow[phiHead[B->omegaIdx[j]]];
			B->psi->cnt++;
		}
	}

#if defined(STOCH_CHECK)
	printf("Deterministic component of reduced cost  = ");
	printSparseVector(B->gBar+1, basisHead, numRows);
#endif

	mem_free(basisHead); mem_free(phiHead); mem_free(basicCost);
	return B;
}//END newBasis()

/* This function allocates a new basisType data structure which holds all the unique basis discovered by the algorithm. It returns a pointer to the
 * structure. */
basisType *newBasisType(string senx, int numIter, int numCols, int numRows, int wordLength) {
	basisType *basis;
	int n;

	if ( !(basis = (basisType *) mem_malloc(sizeof(basisType))))
		errMsg("allocation", "newBasisType", "basis", 0);
	if ( !(basis->vals = (oneBasis **) arr_alloc(numIter, oneBasis *)))
		errMsg("allocation", "newBasisType", "basis->vals", 0);
	if ( !(basis->obsFeasible = (BOOL **) arr_alloc(numIter, BOOL *)))
		errMsg("allocation", "newBasisType", "basis->obsFeasible", 0);
	if ( !(basis->feasSenx = (vector) arr_alloc(numRows+1, double)) )
		errMsg("allocation", "newBasisType", "basis->feasSenx", 0);
	basis->cnt = 0;
	basis->cCodeLen = ceil(numCols/wordLength) + 1;
	basis->rCodeLen = ceil(numRows/wordLength) + 1;

	for (n = 1; n <= numRows; n++ ) {
		if ( senx[n-1] == 'E' )
			basis->feasSenx[n] = 0;
		else if ( senx[n-1] == 'L' )
			basis->feasSenx[n] = -1;
		else if ( senx[n-1] == 'G' )
			basis->feasSenx[n] = 1;
	}

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
		if (B->gBar) mem_free(B->gBar);
		if ( B->phi) {
			for ( n = 0; n < B->phiLength; n++ )
				if (B->phi[n]) mem_free(B->phi[n]);
			mem_free(B->phi);
		}
		if (B->psi) freeSparseMatrix(B->psi);
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
