/*
 * rc.c
 *
 *  Created on: Aug 15, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"

extern configType 	config;
extern ENVptr		env;

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
	if ( newOmegaFlag ) {
		calcDeltaCol(prob->num, prob->coord, cell->lambda, cell->omega->vals[omegaIdx], omegaIdx, cell->delta);
		checkBasisFeasibility(prob->num, prob->coord, cell->basis, cell->omega, newOmegaFlag, omegaIdx, config.MAX_ITER);
	}

	if ( prob->num->rvdOmCnt > 0 ) {
		/* The random variables corresponding to cost coefficients are listed at the end of vector, the offset is used to index them */
		offset = prob->num->rvbOmCnt + prob->num->rvCOmCnt;

		/* If the cost-coefficients are random, update the basis structure. */
		basisIdx = calcBasis(cell->subprob->lp, cell->basis, prob->dBar, cstat, prob->num->cols, rstat, prob->num->rows,
				prob->coord->rvCols, prob->num->rvdOmCnt, &newBasisFlag);

		if ( cell->basis->vals[basisIdx]->phiLength > 0) {
			/* Decompose the dual solution into deterministic and stochastic components. */
			decomposeDualSolution(cell->basis->vals[basisIdx]->phi, cell->omega->vals[omegaIdx]+offset+1, cell->piS,
					cell->basis->vals[basisIdx]->omegaIdx, cell->basis->vals[basisIdx]->phiLength, prob->num->rows);
		}

		if ( newBasisFlag ) {
			/* Establish the feasibility of the new basis with respect to all the observations encountered thus far */
			checkBasisFeasibility(prob->num, prob->coord, cell->basis, cell->omega, FALSE, basisIdx, config.MAX_ITER);

			/* Calculations with respect to deterministic component of the dual solution */
			/* Extract the deterministic component of dual solutions corresponding to rows with random elements in them */
			lambdaIdx = cell->basis->vals[basisIdx]->lambdaIdx[0] = calcLambda(prob->num, prob->coord, cell->piS, cell->lambda, &newLambdaFlag);

			/* Compute the product of deterministic component of dual solution with deterministic (mean value) right-hand side and transfer matrix. */
			cell->basis->vals[basisIdx]->sigmaIdx[0] = calcSigma(prob->num, prob->coord, prob->bBar, prob->Cbar, cell->piS, cell->mubBar,
					lambdaIdx, newLambdaFlag, cell->k, cell->sigma, &newSigmaFlag);

			/* Only need to calculate row if a distinct lambda was found. We could use Pi, instead of lambda(Pi), for this calculation, */
			/* and save the time for expanding/reducing vector even though the lambda is the same, the current Pi might be a
	     distinct one due to the variations in sigma*/
			if (newLambdaFlag)
				calcDeltaRow(config.MAX_ITER, prob->num, prob->coord, cell->omega, cell->lambda, lambdaIdx, cell->delta);

			/* Calculations with respect to stochastic component of the dual solution */
			for (cnt = 0; cnt < cell->basis->vals[basisIdx]->phiLength; cnt++ ) {
				/* Extract the deterministic component of dual solutions corresponding to rows with random elements in them */
				lambdaIdx = cell->basis->vals[basisIdx]->lambdaIdx[cnt+1] = calcLambda(prob->num, prob->coord, cell->basis->vals[basisIdx]->phi[cnt], cell->lambda, &newLambdaFlag);

				/* Compute the product of deterministic component of dual solution with deterministic (mean value) right-hand side and transfer matrix. */
				cell->basis->vals[basisIdx]->sigmaIdx[cnt+1] = calcSigma(prob->num, prob->coord, prob->bBar, prob->Cbar, cell->basis->vals[basisIdx]->phi[cnt], cell->mubBar,
						lambdaIdx, newLambdaFlag, cell->k, cell->sigma, &newSigmaFlag);

				/* Only need to calculate row if a distinct lambda was found. We could use Pi, instead of lambda(Pi), for this calculation, */
				/* and save the time for expanding/reducing vector even though the lambda is the same, the current Pi might be a
		     distinct one due to the variations in sigma*/
				if (newLambdaFlag)
					calcDeltaRow(config.MAX_ITER, prob->num, prob->coord, cell->omega, cell->lambda, lambdaIdx, cell->delta);
			}
		}
	}
	else {
		/* extract the dual solutions corresponding to rows with random elements in them */
		lambdaIdx = calcLambda(prob->num, prob->coord, cell->piS, cell->lambda, &newLambdaFlag);

		/* compute Pi x bBar and Pi x Cbar */
		sigmaIdx = calcSigma(prob->num, prob->coord, prob->bBar, prob->Cbar, cell->piS, cell->mubBar, lambdaIdx, newLambdaFlag, cell->k, cell->sigma, &newSigmaFlag);

		if ( newSigmaFlag ) {
			basisIdx = cell->basis->cnt++;
			cell->basis->vals[basisIdx] = newBasis(0, NULL, NULL, prob->num->cols);
			checkBasisFeasibility(prob->num, prob->coord, cell->basis, cell->omega, FALSE, basisIdx, config.MAX_ITER);

			cell->basis->vals[basisIdx]->lambdaIdx[0] = lambdaIdx;
			cell->basis->vals[basisIdx]->sigmaIdx[0]  = sigmaIdx;
		}
		else
			basisIdx = sigmaIdx;

		/* Only need to calculate row if a distinct lambda was found. We could use Pi, instead of lambda(Pi), for this calculation, */
		/* and save the time for expanding/reducing vector even though the lambda is the same, the current Pi might be a
			     distinct one due to the variations in sigma*/
		if (newLambdaFlag)
			calcDeltaRow(config.MAX_ITER, prob->num, prob->coord, cell->omega, cell->lambda, lambdaIdx, cell->delta);
	}

	mem_free(cstat);
	mem_free(rstat);
	return basisIdx;
}//End stochasticUpdates()

int calcBasis(LPptr lp, basisType *basis, sparseVector *dBar, intvec cstat, int numCols, intvec rstat, int numRows, intvec rvCols, int rvdOmCnt, BOOL *newBasisFlag) {
	unsigned long *codedCol, *codedRow;
	vector	costVector, tempPsiRow;
	intvec 	basisHead, randBasisHead;
	int		cnt, i, j;

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

	/* Allocate memory for the basis header. */
	if ( !(basisHead = (intvec) arr_alloc(numRows+1, int)) )
		errMsg("allocation", "calcBasis", "basisHead", 0);
	if ( !(randBasisHead = (intvec) arr_alloc(numRows+1, int)) )
		errMsg("allocation", "calcBasis", "randBasisHead", 0);

	/* New basis encountered, add it to the list */
	(*newBasisFlag) = TRUE;
	basis->vals[cnt] = newBasis(numRows, codedCol, codedRow, numCols);

	/* Compute the phi matrix associated with the current basis. We begin by first identifying the basis header. A negative value in basis header indicates a slack row. */
	getBasisHead(lp, basisHead+1, NULL);

	/* Compute the phi matrix header and extract the basis (of the dual) inverse matrix rows corresponding to the header. */
	for ( i = 1; i <= numRows; i++ ) {
		if ( basisHead[i] >= 0 ) {
			/* corresponds to a basic row */
			j = 1;
			while ( j <= rvdOmCnt ) {
				if ( (rvCols[j] - 1) == basisHead[i] ) {
					/* basis column with random cost coefficient */
					randBasisHead[basis->vals[cnt]->phiLength+1] = i;
					basis->vals[cnt]->phiHeader[basis->vals[cnt]->phiLength+1] = basisHead[i]+1;
					basis->vals[cnt]->omegaIdx[basis->vals[cnt]->phiLength+1] = j;
					if ( !(basis->vals[cnt]->phi[basis->vals[cnt]->phiLength] = (vector) arr_alloc(numRows+1, double)) )
						errMsg("allocation", "calcBasis", "basis->vals[cnt]->phi[i]", 0);
					getBasisInvRow(lp, i-1, basis->vals[cnt]->phi[basis->vals[cnt]->phiLength]+1);
					basis->vals[cnt]->phiLength++;
				}
				j++;
			}
		}
	}

	if ( basis->vals[cnt]->phiLength > 0 ) {
		/* reallocate memory to phi header and matrix */
		basis->vals[cnt]->phiHeader = (intvec) mem_realloc(basis->vals[cnt]->phiHeader, (basis->vals[cnt]->phiLength+1)*sizeof(int));
		basis->vals[cnt]->phi 		= (vector *) mem_realloc(basis->vals[cnt]->phi, basis->vals[cnt]->phiLength*sizeof(vector));
		basis->vals[cnt]->omegaIdx = (intvec) mem_realloc(basis->vals[cnt]->omegaIdx, (basis->vals[cnt]->phiLength+1)*sizeof(int));
		basis->vals[cnt]->lambdaIdx = (intvec) mem_realloc(basis->vals[cnt]->lambdaIdx, (basis->vals[cnt]->phiLength+1)*sizeof(int));
		basis->vals[cnt]->sigmaIdx = (intvec) mem_realloc(basis->vals[cnt]->sigmaIdx, (basis->vals[cnt]->phiLength+1)*sizeof(int));

		/* Initialize elements which will be used for feasibility check */
		if ( !(basis->vals[cnt]->psi = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))) )
			errMsg("allocation", "newBasis", "B->psi", 0);
		basis->vals[cnt]->psi->val = (vector) arr_alloc(numCols*basis->vals[cnt]->phiLength+1, double);
		basis->vals[cnt]->psi->col = (intvec) arr_alloc(numCols*basis->vals[cnt]->phiLength+1, int);
		basis->vals[cnt]->psi->row = (intvec) arr_alloc(numCols*basis->vals[cnt]->phiLength+1, int);
		basis->vals[cnt]->psi->cnt = 0;
	}
	else {
		/* If the basis does not include columns that do not have random cost coefficients. */
		mem_free(basis->vals[cnt]->phiHeader); 	basis->vals[cnt]->phiHeader = NULL;
		mem_free(basis->vals[cnt]->phi);		basis->vals[cnt]->phi 		= NULL;
		mem_free(basis->vals[cnt]->omegaIdx);	basis->vals[cnt]->omegaIdx  = NULL;
		basis->vals[cnt]->lambdaIdx = (intvec) mem_realloc(basis->vals[cnt]->lambdaIdx, sizeof(int));
		basis->vals[cnt]->sigmaIdx = (intvec) mem_realloc(basis->vals[cnt]->sigmaIdx, sizeof(int));

		/* Feasibility in this case is assessed using only the _g_ vector in basis structure */
		basis->vals[cnt]->psi = NULL;
	}


	/* Extract the basic variable cost vector and psi matrix (the tableau entries) */
	costVector = expandVector(dBar->val, dBar->col, dBar->cnt, numCols);
	if ( !(tempPsiRow = (vector) arr_alloc(numRows+1, double)) )
		errMsg("allocation", "calcBasis", "tempPsiRow", 0);

	for ( i = 1; i <= numCols; i++ ) {
		getBasisInvACol(lp, i-1, tempPsiRow+1);

		basis->vals[cnt]->g[i] = costVector[i];
		for ( j = 1; j <= numRows; j++ )
			basis->vals[cnt]->g[i] -= tempPsiRow[j]*costVector[basisHead[j]+1];

		for ( j = 1; j <= basis->vals[cnt]->phiLength; j++ ) {
			basis->vals[cnt]->psi->val[basis->vals[cnt]->psi->cnt+1]   = tempPsiRow[randBasisHead[j]];
			basis->vals[cnt]->psi->row[basis->vals[cnt]->psi->cnt+1]   = i;
			basis->vals[cnt]->psi->col[basis->vals[cnt]->psi->cnt+1] = basis->vals[cnt]->phiHeader[j];
			basis->vals[cnt]->psi->cnt++;
		}
	}

#if defined (STOCH_CHECK)
	printf("New basis identified     :: %d\n", cnt);
	printf("\tNumber of basic stochastic columns = %d\n", basis->vals[cnt]->phiLength);
	if ( basis->vals[cnt]->phiLength > 0 ) {
		printf("\tBasic stochastic columns           = "); printIntvec(basis->vals[cnt]->phiHeader, basis->vals[cnt]->phiLength, NULL);
		printf("\tStochastic cost variable           = "); printIntvec(basis->vals[cnt]->omegaIdx, basis->vals[cnt]->phiLength, NULL);
		printf("\tPhi = ");
		for (i = 0; i < basis->vals[cnt]->phiLength; i++ ) {
			printf("\t\t"); printVector(basis->vals[cnt]->phi[i], numRows, NULL);
		}
	}
	else {
		printf("\tBasic stochastic columns           = NULL\n");
		printf("\tStochastic cost variable           = NULL\n");
		printf("\tPhi                                = NULL\n");
	}
#endif

	mem_free(costVector); mem_free(basisHead); mem_free(tempPsiRow);
	return basis->cnt++;

}//END calcBasis()

int checkBasisFeasibility(numType *num, coordType *coord, basisType *basis, omegaType *omega, BOOL newOmegaFlag, int elemIdx, int maxIter) {
	int n, c, offset;
	vector costVector;
	sparseVector cOmega;

	offset = num->rvbOmCnt + num->rvCOmCnt;
	cOmega.cnt = num->rvdOmCnt; cOmega.col = coord->omegaCol+offset;

	if ( !(costVector = (vector) arr_alloc(num->cols+1, double)) )
		errMsg("allocation", "calcBasisFeasibility", "costVector", 0);

	if ( newOmegaFlag ) {
		cOmega.val = omega->vals[elemIdx] + offset;
		for ( n = 0; n < basis->cnt; n++ ) {
			basis->obsFeasible[n][elemIdx] = TRUE;
			if  ( cOmega.cnt > 0 ) {
				copyVector(basis->vals[n]->g, costVector, num->cols, TRUE);
				addVectors(costVector, cOmega.val, cOmega.col, cOmega.cnt);
				if ( basis->vals[n]->phiLength > 0 ) {
					MSparsexvSub(basis->vals[n]->psi, cOmega.val, costVector);
				}
				c = 1;
				while ( c <= num->cols) {
					if ( costVector[c] < 0 )
						basis->obsFeasible[n][elemIdx] = FALSE;
					c++;
				}
			}
		}
	}
	else {
		if ( !(basis->obsFeasible[elemIdx] = (BOOL*) arr_alloc(maxIter, BOOL)) )
			errMsg("allocation", "calcBasisFeasibility", "basis->obsFeasibility[n]", 0);
		for ( n = 0; n < omega->cnt; n++ ) {
			basis->obsFeasible[elemIdx][n] = TRUE;
			if ( cOmega.cnt > 0 ) {
				cOmega.val = omega->vals[n] + offset;
				copyVector(basis->vals[elemIdx]->g, costVector, num->cols, TRUE);
				addVectors(costVector, cOmega.val, cOmega.col, cOmega.cnt);
				if ( basis->vals[elemIdx]->phiLength > 0 ) {
					MSparsexvSub(basis->vals[elemIdx]->psi, cOmega.val, costVector);
				}
				c = 1;
				while ( c <= num->cols) {
					if ( costVector[c] < 0 )
						basis->obsFeasible[elemIdx][n] = FALSE;
					c++;
				}
			}
		}
	}

	mem_free(costVector);
	return 0;
}//END calcBasisFeasibility()

int decomposeDualSolution(vector *phi, vector omegaVals, vector Pi, intvec phiOmegaIdx, int phiLength, int numRows) {
	int n, i;

	for ( n = 0; n < phiLength; n++ )
		for ( i = 1; i <= numRows; i++ )
			Pi[i] -= phi[n][i]*omegaVals[n];

	return 0;
}//END decomposeDualSolution()

/* The function encodes an integer vector _stream_ of given length _len_ into an unsigned long vector _codeWord_ */
unsigned long *encodeIntvec(intvec stream, int len, int wordLength) { /* TODO: After merging 2SD_randomCost branch into main, move this subroutine to utilities */
	unsigned long *codeWord, temp;
	int j, group, shift, codeLength;

	codeLength = ceil((double) len/ (double) wordLength) + 1;

	if ( !(codeWord = (unsigned long *) arr_alloc(codeLength, unsigned long)))
		errMsg("allocation", "encodeIntVec", "codeWord", 0);

	for (j = 1; j <= len; j++) {
		group = j/wordLength + 1;
		shift = wordLength - (j % wordLength);
		temp = (unsigned long) stream[j] << shift;
		codeWord[group] |= temp;
	}

	return codeWord;
}//END encodeIntVec()

BOOL equalLongIntvec(unsigned long *a, unsigned long *b, int len) {  /* TODO: After merging 2SD_randomCost branch into main, move this subroutine to utilities */
	int		cnt;

	for (cnt = 1; cnt <= len; cnt++)
		if ( a[cnt] != b[cnt] )
			return FALSE;

	return TRUE;
}//END equalLongIntvec()

int getBasisHead(LPptr lp, intvec head, vector basicX) {
	int status;

	status = CPXgetbhead(env, lp, head, basicX);
	if ( status )
		solverErrmsg(status);

	return status;
}//END getBasisHead()

int getBasisInvRow(LPptr lp, int i, vector phi) {
	int status;

	status = CPXbinvrow(env, lp, i, phi);
	if ( status )
		solverErrmsg(status);

	return status;
}//END getBasicInvRow()

int getBasisInvCol(LPptr lp, int i, vector phi) {
	int status;

	status = CPXbinvcol(env, lp, i, phi);
	if ( status )
		solverErrmsg(status);

	return status;
}//END getBasicInvRow()

int getBasisInvARow(LPptr lp, int i, vector phi) {
	int status;

	status = CPXbinvarow(env, lp, i, phi);
	if ( status )
		solverErrmsg(status);

	return status;
}//END getBasicInvRow()

int getBasisInvACol(LPptr lp, int i, vector phi) {
	int status;

	status = CPXbinvacol(env, lp, i, phi);
	if ( status )
		solverErrmsg(status);

	return status;
}//END getBasicInvRow()

/* This subroutine extracts elements which are common to the two input integer vectors _a_ and _b_ */
intvec intvecIntersect(intvec a, intvec b, int lenA, int lenB) {
	intvec inter;
	int	cnt, n;

	if ( !(inter = (intvec) arr_alloc(max(lenA, lenB)+1, int)) )
		errMsg("allocation", "intvecIntersect", "inter", 0);

	cnt = 1;
	for ( n = 1; n <= lenA; n++ )
		if ( isElementIntvec(b, lenB, a[n]) )
			inter[cnt++] = a[n];

	return inter;

}//END intvecIntersect()

/* This subroutine checks to see if a integer scalar is an element of integer vector. If so, the subroutine will return the index. If not, a value of -1 is returned. */
int isElementIntvec(intvec vec, int lenVec, int elem) {
	int n = 1;

	while ( vec[n] != elem && n <= lenVec )
		n++;

	if ( n == (lenVec+1) )
		return -1;
	else
		return n;

}

void subVectors(vector a, vector b, intvec indices, int len){
	int n;

	if ( indices == NULL ) {
		for ( n = 1; n <= len; n++ )
			a[n] -= b[n];
	}
	else {
		for ( n = 1; n <= len; n++ )
			a[indices[n]] -= b[n];
	}
	a[0] = oneNorm(a+1, len);

}//END copy_arr()
