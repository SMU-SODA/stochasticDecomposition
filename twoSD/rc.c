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
		return 1;
	}
	/* Record the dual and reduced cost on bounds. */
	if ( getDual(cell->subprob->lp, cell->piS, prob->num->rows) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to get the dual", 0);
		return 1;
	}
	if ( computeMU(cell->subprob->lp, cstat,  prob->num->cols, &cell->mubBar) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to compute mubBar for subproblem", 0);
		return 1;
	}

	/* Update the column of delta structure if a new observation was encountered. */
	if ( newOmegaFlag )
		calcDeltaCol(prob->num, prob->coord, cell->lambda, cell->omega->vals[omegaIdx], omegaIdx, cell->delta);

	if ( prob->num->rvdOmCnt > 0 ) {
		/* The random variables corresponding to cost coefficients are listed at the end of vector, the offset is used to index them */
		offset = prob->num->rvbOmCnt + prob->num->rvCOmCnt;

		/* If the cost-coefficients are random, update the basis structure. */
		basisIdx = calcBasis(cell->basis, cell->subprob->lp, cstat, prob->num->cols, rstat, prob->num->rows,
				prob->coord->rvCols+offset, prob->num->rvdOmCnt, &newBasisFlag);

		if ( cell->basis->vals[basisIdx]->phiLength > 0) {
			/* Decompose the dual solution into deterministic and stochastic components. */
			decomposeDualSolution(cell->basis->vals[basisIdx]->phi, cell->omega->vals[omegaIdx]+offset, cell->piS,
					cell->basis->vals[basisIdx]->omegaIdx, cell->basis->vals[basisIdx]->phiLength, prob->num->rows);
		}

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
	else {
		/* extract the dual solutions corresponding to rows with random elements in them */
		lambdaIdx = calcLambda(prob->num, prob->coord, cell->piS, cell->lambda, &newLambdaFlag);

		/* compute Pi x bBar and Pi x Cbar */
		sigmaIdx = calcSigma(prob->num, prob->coord, prob->bBar, prob->Cbar, cell->piS, cell->mubBar, lambdaIdx, newLambdaFlag, cell->k, cell->sigma, &newSigmaFlag);

		if ( newSigmaFlag ) {
			cell->basis->vals[cell->basis->cnt] 			  = newBasis(0, NULL, NULL);
			cell->basis->vals[cell->basis->cnt]->lambdaIdx[0] = lambdaIdx;
			cell->basis->vals[cell->basis->cnt]->sigmaIdx[0]  = sigmaIdx;
		}

		/* Only need to calculate row if a distinct lambda was found. We could use Pi, instead of lambda(Pi), for this calculation, */
		/* and save the time for expanding/reducing vector even though the lambda is the same, the current Pi might be a
			     distinct one due to the variations in sigma*/
		if (newLambdaFlag)
			calcDeltaRow(config.MAX_ITER, prob->num, prob->coord, cell->omega, cell->lambda, lambdaIdx, cell->delta);
	}

	mem_free(cstat); mem_free(rstat);
	return 0;
}//End stochasticUpdates()

int calcBasis(basisType *basis, LPptr lp, intvec cstat, int numCols, intvec rstat, int numRows, intvec rvCols, int rvdOmCnt, BOOL *newBasisFlag) {
	unsigned long *codedCol, *codedRow;
	intvec 	basisHead;
	int		cnt, i, j;

	/* allocate memory */
	if ( !(basisHead = (intvec) arr_alloc(numRows+1, int)) )
		errMsg("allocation", "calcBasis", "basisHead", 0);

	/* encode the row and column status */
	codedCol = encodeIntvec(cstat, numCols, WORDLENGTH);
	codedRow = encodeIntvec(rstat, numRows, WORDLENGTH);

	/* check to see if the current basis was encountered before */
	for ( cnt = 0; cnt < basis->cnt; cnt++ ) {
		if ( equalLongIntvec(codedCol, basis->vals[cnt]->cCode, basis->cCodeLen) && equalLongIntvec(codedRow, basis->vals[cnt]->rCode, basis->rCodeLen) ) {
			/* The basis is the same as one encountered before */
			basis->vals[cnt]->weight++;
			mem_free(codedRow); mem_free(codedCol);
			(*newBasisFlag) = FALSE;
			return cnt;
		}
	}

	/* New basis encountered, add it to the list */
	(*newBasisFlag) = TRUE;
	basis->vals[cnt] = newBasis(numRows, codedCol, codedRow);

	/* Compute the phi matrix associated with the current basis. We begin by first identifying the basis header. A negative value in basis header indicates a slack row. */
	getBasisHead(lp, basisHead+1, NULL);

	/* Compute the phi matrix header and extract the basis inverse matrix rows corresponding to the header */
	for ( i = 1; i <= numRows; i++ ) {
		if ( basisHead[i] >= 0 ) {
			/* corresponds to a basic row */
			j = 1;
			while ( j <= rvdOmCnt ) {
				if ( (rvCols[j] - 1) == basisHead[i] ) {
					/* basis column with random cost coefficient */
					basis->vals[cnt]->phiHeader[basis->vals[cnt]->phiLength] = basisHead[i];
					basis->vals[cnt]->omegaIdx[basis->vals[cnt]->phiLength] = j;
					if ( !(basis->vals[cnt]->phi[basis->vals[cnt]->phiLength] = (vector) arr_alloc(numRows+1, double)) )
						errMsg("allocation", "calcBasis", "basis->vals[cnt]->phi[i]", 0);
					getBasisInvRow(lp, i, basis->vals[cnt]->phi[basis->vals[cnt]->phiLength]+1);
					basis->vals[cnt]->phiLength++;
				}
				j++;
			}
		}
	}

	if ( basis->vals[cnt]->phiLength > 0 ) {
		/* reallocate memory to phi header and matrix */
		basis->vals[cnt]->phiHeader = (intvec) mem_realloc(basis->vals[cnt]->phiHeader, basis->vals[cnt]->phiLength*sizeof(int));
		basis->vals[cnt]->phi 		= (vector *) mem_realloc(basis->vals[cnt]->phi, basis->vals[cnt]->phiLength*sizeof(vector));
		basis->vals[cnt]->omegaIdx = (intvec) mem_realloc(basis->vals[cnt]->omegaIdx, basis->vals[cnt]->phiLength*sizeof(int));
		basis->vals[cnt]->lambdaIdx = (intvec) mem_realloc(basis->vals[cnt]->lambdaIdx, (basis->vals[cnt]->phiLength+1)*sizeof(int));
		basis->vals[cnt]->sigmaIdx = (intvec) mem_realloc(basis->vals[cnt]->sigmaIdx, (basis->vals[cnt]->phiLength+1)*sizeof(int));
	}
	else {
		/* If the basis does not include columns that do not have random cost coefficients. */
		mem_free(basis->vals[cnt]->phiHeader); 	basis->vals[cnt]->phiHeader = NULL;
		mem_free(basis->vals[cnt]->phi);		basis->vals[cnt]->phi 		= NULL;
		mem_free(basis->vals[cnt]->omegaIdx);	basis->vals[cnt]->omegaIdx  = NULL;
		basis->vals[cnt]->lambdaIdx = (intvec) mem_realloc(basis->vals[cnt]->lambdaIdx, sizeof(int));
		basis->vals[cnt]->sigmaIdx = (intvec) mem_realloc(basis->vals[cnt]->sigmaIdx, sizeof(int));
	}

	mem_free(basisHead);
	return basis->cnt++;

}//END calcBasis()

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

	codeLength = ceil(len/wordLength) + 1;

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
