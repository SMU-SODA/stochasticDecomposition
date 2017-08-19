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

extern ENVptr	env;

int calcBasis(LPptr lp, int numCols, int numRows, int numCostRVs, intvec costRVcols, basisType *basis, BOOL *newBasisFlag) {
	intvec cstat, rstat, basisHead;;
	unsigned long *codedRow, *codedCol;
	int    i, j, cnt;

	/* allocate memory */
	if ( !(cstat = (intvec) arr_alloc(numCols+1, int)))
		errMsg("allocation", "getIndexNumber", "cstat", 0);
	if ( !(rstat = (intvec) arr_alloc(numRows+1, int)))
		errMsg("allocation", "getIndexNumber", "rstat", 0);
	if ( !(basisHead = (intvec) arr_alloc(numRows+1, int)) )
		errMsg("allocation", "calcBasis", "basisHead", 0);


	/* obtain the status of columns and rows in the basis */
	getBasis(lp, cstat+1, rstat+1);  /* TODO: can incorporate computeMU into this program. */

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
	basis->vals[cnt] = newBasis(numRows);
	basis->vals[cnt]->cCode  = codedCol;
	basis->vals[cnt]->rCode  = codedRow;
	basis->vals[cnt]->weight = 1;

	/* Compute the phi matrix associated with the current basis. We begin by first identifying the basis header. A negative value in basis header indicates a slack row. */
	getBasisHead(lp, basisHead+1, NULL);

	/* Compute the phi matrix header and extract the basis inverse matrix rows corresponding to the header */
	for ( i = 1; i <= numRows; i++ ) {
		if ( basisHead[i] >= 0 ) {
			/* corresponds to a basic row */
			j = 1;
			while ( j <= numCostRVs ) {
				if ( (costRVcols[j] - 1) == basisHead[i] ) {
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

	/* reallocate memory to phi header and matrix */
	if ( basis->vals[cnt]->phiLength ) {
		basis->vals[cnt]->phiHeader = (intvec) mem_realloc(basis->vals[cnt]->phiHeader, basis->vals[cnt]->phiLength*sizeof(int));
		basis->vals[cnt]->phi 		= (vector *) mem_realloc(basis->vals[cnt]->phi, basis->vals[cnt]->phiLength*sizeof(vector));
	}
	else {
		/* If the basis does not include columns that do not have random cost coefficients. */
		basis->vals[cnt]->phiHeader = NULL;
		basis->vals[cnt]->phi 		= NULL;
	}

	mem_free(basisHead);
	mem_free(cstat); mem_free(rstat);

	return basis->cnt++;
}//End calcBasis()

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
