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
	int 	basisIdx;
	BOOL	newBasisFlag;

	/* Allocate memory. */
	if ( !(cstat = (intvec) arr_alloc( prob->num->cols+1, int)))
		errMsg("allocation", "getIndexNumber", "cstat", 0);
	if ( !(rstat = (intvec) arr_alloc( prob->num->rows+1, int)))
		errMsg("allocation", "getIndexNumber", "rstat", 0);

	/* Obtain the status of columns and rows in the basis. */
	if ( getBasis(cell->subprob->lp, cstat+1, rstat+1) ) {
		errMsg("algorithm", "stochastics", "failed to get the basis column and row status", 0);
		return 1;
	}
	/* Record the dual and reduced cost on bounds. */
	if ( getDual(cell->subprob->lp, cell->piS, prob->num->rows) ) {
		errMsg("algorithm", "stochastics", "failed to get the dual", 0);
		return 1;
	}
	if ( computeMU(cell->subprob->lp, cstat,  prob->num->cols, &cell->mubBar) ) {
		errMsg("algorithm", "stochastics", "failed to compute mubBar for subproblem", 0);
		return 1;
	}

	/* Update the column of delta structure if a new observation was encountered. */
	if ( newOmegaFlag )
		calcDeltaCol(prob->num, prob->coord, cell->lambda, cell->omega->vals[omegaIdx], omegaIdx, cell->delta);

	if ( prob->num->rvdOmCnt > 0 ) {
		/* If the cost-coefficients are random, update the basis structure. */
		basisIdx = calcBasis(cell->basis, cell->subprob->lp, cstat, prob->num->cols, rstat, prob->num->rows, prob->coord->rvCols, prob->num->rvdOmCnt, &newBasisFlag);

		if ( cell->basis->vals[basisIdx]->phiLength > 0) {
			/* TODO: Decompose the dual solution into deterministic and stochastic components. */

		}

		// TODO: If there is a new basis, then the following need to be computed. */
		calcLambdaSigma(prob->num, prob->coord, prob->bBar, prob->Cbar, cell->piS, cell->mubBar, cell->k, cell->basis, basisIdx, &newBasisFlag);
	}
	else {
		/* The cost coefficients do not have randomness. Update lambda and sigma as based on the dual solutions alone */
		calcLambdaSigma(prob->num, prob->coord, prob->bBar, prob->Cbar, cell->piS, cell->mubBar, cell->k, cell->basis, -1, &newBasisFlag);
	}

	/* Only need to calculate row if a distinct lambda was found. We could use Pi, instead of lambda(Pi), for this calculation, */
	/* and save the time for expanding/reducing vector even though the lambda is the same, the current Pi might be a
     distinct one due to the variations in sigma*/
	if (newBasisFlag)
		calcDeltaRow(config.MAX_ITER, prob->num, prob->coord, cell->omega, cell->lambda, basisIdx, cell->delta);


	mem_free(cstat); mem_free(rstat);
	return 0;
}//End stochastics()

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
					if ( !(basis->vals[cnt]->phi[basis->vals[cnt]->phiLength] = (vector) arr_alloc(numRows, double)) )
						errMsg("allocation", "calcBasis", "basis->vals[cnt]->phi[i]", 0);
					getBasisInvRow(lp, i, basis->vals[cnt]->phi[basis->vals[cnt]->phiLength]);
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
//		basis->vals[cnt]
	}
	else {
		/* If the basis does not include columns that do not have random cost coefficients. */
		mem_free(basis->vals[cnt]->phiHeader); 	basis->vals[cnt]->phiHeader = NULL;
		mem_free(basis->vals[cnt]->phi);		basis->vals[cnt]->phi 		= NULL;
	}

	mem_free(basisHead);
	return basis->cnt++;

}//END calcBasis()

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
