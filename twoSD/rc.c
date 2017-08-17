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

int calcBasis(LPptr lp, int numCols, int numRows, basisType *basis, BOOL *newBasisFlag) {
	intvec cstat, rstat;
	unsigned long *codedRow, *codedCol;
	int    cnt;

	/* allocate memory */
	if ( !(cstat = (intvec) arr_alloc(numCols+1, int)))
		errMsg("allocation", "getIndexNumber", "cstat", 0);
	if ( !(rstat = (intvec) arr_alloc(numRows+1, int)))
		errMsg("allocation", "getIndexNumber", "rstat", 0);

	/* obtain the status of columns and rows in the basis */
	getBasis(lp, cstat+1, rstat+1);  /* TODO: can incorporate computeMU into this program. */

	/* encode the row and column status */
	codedCol = encodeIntvec(cstat, numCols, WORDLENGTH);
	codedRow = encodeIntvec(rstat, numRows, WORDLENGTH);

	/* check to see if the current basis was encountered before */
	for ( cnt = 0; cnt < basis->cnt; cnt++ ) {
		if ( equalLongIntvec(codedCol, basis->cCode[cnt], basis->cCodeLen) && equalLongIntvec(codedRow, basis->rCode[cnt], basis->rCodeLen) ) {
			/* The basis is the same as one encountered before */
			basis->weight[cnt]++;
			mem_free(codedRow); mem_free(codedCol);
			(*newBasisFlag) = FALSE;
			return cnt;
		}
	}

	/* New basis encountered, add it to the list */
	(*newBasisFlag) = TRUE;
	basis->cCode[cnt] = codedCol;
	basis->rCode[cnt] = codedRow;
	basis->weight[cnt] = 1;

	mem_free(cstat); mem_free(rstat);

	return basis->cnt++;
}//End getIndexNumber()

/* The function encodes an integer vector _stream_ of given length _len_ into an unsigned long vector _codeWord_ */
unsigned long *encodeIntvec(intvec stream, int len, int wordLength) { /* TODO: After merging 2SD_randomCost branch into main, move this subroutine to utilities */
	unsigned long *codeWord, temp;
    int j, group, shift;

    if ( !(codeWord = (unsigned long *) arr_alloc(len/wordLength + 2, unsigned long)))
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
