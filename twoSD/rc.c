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

extern ENVptr		env;

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
