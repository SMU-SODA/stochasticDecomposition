/*
 * utility.c
 *
 *  Created on: Apr 20, 2014
 *      Author: Harsha Gangammanavar
 */

#include "utils.h"

extern long MEM_USED;

FILE *openFile(cString outputDir, cString fname, char *mode) {
	FILE *fptr;
	char buffer[2*BLOCKSIZE];

	strcpy(buffer, outputDir);
	strcat(buffer, fname);

	fptr = fopen(buffer, mode);
	if ( fptr == NULL ) {
		fprintf(stderr, "failed to open file %s", fname);
		return NULL;
	}
	return fptr;
}//END openFile()

void createOutputDir(cString outputDir, cString algoName, cString probName) {
	struct stat st;
	char buffer[2*BLOCKSIZE];

	strcat(outputDir,algoName);
	strcat(outputDir,"/");
	if ( stat(outputDir, &st) ) {
		sprintf(buffer, "mkdir %s", outputDir);
		system(buffer); 	fflush(stdout);
	}
	strcat(outputDir, probName);
	strcat(outputDir, "/");
	if ( stat(outputDir, &st) ) {
		sprintf(buffer, "mkdir %s", outputDir);
		system(buffer); 	fflush(stdout);
	}
	else {
		sprintf(buffer, "rm -r %s*", outputDir);
		system(buffer); 	fflush(stdout);
	}

}//END createOutputDir()

void errMsg(cString type, cString place, cString item, int quit){
	fprintf(stderr, "\nError :: Type - %s;  Function - %s(); Item - %s\n", type, place, item);
	if (quit)
		exit(1);
}//err_msg()

/* The function getLine() reads an input line, checks its length, and determines the appropriate parsing function to call. In the case of a line consisting of only a single carriage return, the function
 * continues to read in lines until a cString of non_zero length or an EOF is encountered. Upon return of fields from the parsing functions, each field will be sent to remove_spaces for the removal of all
 * blank spaces.  The getLine function then returns control to the calling function. */
int getLine(FILE **input, cString *fields, char *type, int *numFields) {
	char 	input_str[10*BLOCKSIZE], *strptr, *token;
	long	len = 0, n, stat;

	input_str[0] = '*';

	/*  test for and skip over empty lines and comments  */
	while (len <= 1 || input_str[0] == '*') {
		strptr = fgets(input_str, 1000, *input);
		if (strptr == NULL)
			return 1;
		len = strlen(input_str);
	}

	/*  identify type of input cString (title or field)  */
	if (input_str[0] >= '0' && input_str[0] <= 'Z') {
		/* input cString is title cString, marked by type = 't' */
		sscanf(input_str, "%s %s", fields[0], fields[1]);
		type[0] = 't';
		n = 2;
	}
	else {
		/* input cString is field cString, marked by type = 'f' */
		token = strtok(input_str, " \t");
		n = 0;
		while( token != NULL ) {
			strcpy(fields[n++], token);
			token = strtok(NULL, " \t");
		}
		type[0] = 'f';
	}

	(*numFields) = n;
	for ( n = 0; n < (*numFields); n++ ) {
		stat = removeSpaces(fields[n]);
		if ( stat == 0 ) {
			(*numFields)--;
		}
	}

	return 0;
}//END getLine()

/* The function removeSpace() removes additional spaces from a	cString. This function will be used after an input cString has been broken into its appropriate fields according to column values. */
int removeSpaces (char *field) {
	char *p, *q, len = 0;;

	p = field;
	q = p;
	while (*p != '\0' && *p != '\n') {
		if (*p >= 33) {
			*q++ = *p++;
			len++;
		}
		else
			p++;
	}
	*q = '\0';

	return len;
}//END removeSpaces

void trPrint(cString routine, int type){
	if ( type == 1 )
		printf("Entering :: %s()\n", routine);
	else
		printf("Exiting  :: %s()\n", routine);
}//END trPrint()

void *log_alloc(char *cString, void *return_ptr, int size) {
	MEM_USED += size;
	return return_ptr;
}//END log_alloc()


void *log_realloc(char *cString, void *free_ptr, void *alloc_ptr, int size) {
	MEM_USED += size;
	return alloc_ptr;
}//END log_realloc()

void mem_free ( void *ptr ) {

	if ( ptr != NULL ) {
		free(ptr);
		ptr = NULL;
	}
}

double str2float(char *cString){
	double val;
	sscanf(cString, "%lf", &val);
	return val;
}//END str_to_float()

int str2int(char *cString) {
	int val;
	sscanf(cString, "%d", &val);
	return val;
}//END str2int()

/* This function returns the minimum number of bits needed to represent a given number. */
int getNumBits(int num) {
	int 	hi_bit = 1;
	int 	numBits;

	for (numBits = 0; hi_bit <= num; numBits++)
		hi_bit = hi_bit << 1;

	return numBits;
}//END getNumBits()


double oneNorm(dVector a, int len) {
	int		cnt;
	double	sum;

	sum = 0.0;
	for (cnt = 0; cnt < len; cnt++)
		sum += DBL_ABS(a[cnt]);

	return sum;
}//END oneNorm()

double twoNorm(dVector a, dVector b, int len) {
	int 	cnt;
	double	norm = 0.0;

	if (b != NULL)
		for (cnt = 1; cnt <= len; cnt++ )
			norm += pow((a[cnt]-b[cnt]), 2);
	else
		for (cnt = 1; cnt <= len; cnt++ )
			norm += pow(a[cnt], 2);

	norm = sqrt(norm);
	return norm;
}//END twoNorm()

double pNorm(dVector a, dVector b, int len, int p) {
	int 	cnt;
	double	norm = 0.0;

	if (b != NULL)
		for (cnt = 1; cnt <= len; cnt++ )
			norm += pow(fabs(a[cnt]-b[cnt]), p);
	else
		for (cnt = 1; cnt <= len; cnt++ )
			norm += abs(pow(a[cnt], p));

	norm = fabs(pow(norm, 1/(double) p));

	return norm;
}//END twoNorm()

double infNorm(dVector a, dVector b, int len) {
	double norm = 0.0;

	if (b != NULL)
		for (int cnt = 1; cnt <= len; cnt++ )
			norm = maximum(fabs(a[cnt]-b[cnt]), norm);
	else
		for (int cnt = 1; cnt <= len; cnt++ )
			norm = maximum(fabs(a[cnt]), norm);

	return norm;
}//END infNorm()

void calcMeanVariance(dVector x, int lenX, double *mean, double *variance) {
	double 	temp;
	int 	cnt;

	temp = 0.0;
	(*variance) = 0.0; (*mean) = x[0];
	for (cnt = 1; cnt < lenX; cnt++) {
		temp = (*mean);
		(*mean) = (*mean) + (x[cnt] - (*mean)) / (double) (cnt + 1);
		(*variance) = (1 - 1 / (double) cnt) * (*variance) + (cnt + 1) * ((*mean) - temp) * ((*mean) - temp);
	}

}//END calcVariance()

double vXv(dVector a, dVector b, iVector idxCol, int len) {
	double ans = 0.0;
	int n;

	if(idxCol == NULL)
		for ( n = 1; n <= len; n++ )
			ans += a[n]*b[n];
	else
		for (n = 1; n <= len; n++ )
			ans += a[n]*b[idxCol[n]];

	return ans;
}//END vXv()

double vXvSparse(dVector v, sparseVector *vSparse){
	int		cnt;
	double 	ans;

	ans = 0.0;
	for (cnt = 1; cnt <= vSparse->cnt; cnt++)
		ans += vSparse->val[cnt] * v[vSparse->col[cnt]];

	return ans;
}//END vXvSparse()

dVector MSparsexvAdd(sparseMatrix *M, dVector v, dVector ans){
	int	cnt;

	for (cnt = 1; cnt <= M->cnt; cnt++)
		ans[M->row[cnt]] += M->val[cnt] * v[M->col[cnt]];

	return ans;
}//END MSparsexv()

dVector MSparsexvSub(sparseMatrix *M, dVector v, dVector ans){
	int	cnt;

	for (cnt = 1; cnt <= M->cnt; cnt++)
		ans[M->row[cnt]] -= M->val[cnt] * v[M->col[cnt]];

	return ans;
}//END MSparsexvSub()

dVector vxMSparse(dVector v, sparseMatrix *M, int len) {
	int		cnt;
	dVector	ans;

	if(!(ans = (dVector) arr_alloc(len+1, double)))
		errMsg("allocation", "vxMSparse", "ans", 1);

	for (cnt = 1; cnt <= M->cnt; cnt++)
		ans[M->col[cnt]] += v[M->row[cnt]] * M->val[cnt];
	ans[0] = oneNorm(ans+1, len);

	return ans;
}//END PIxT()

void vPlusv(dVector a, dVector b, double mult, int len){
	int 	cnt;

	for ( cnt = 1; cnt <= len; cnt++ )
		a[cnt] += mult*b[cnt];
	a[0] = oneNorm(a+1, len);

}//END vPlusv()

double smooth(double new, double old, double factor) {
	return factor*new + (1-factor)*old;
}//END smooth();

dVector reduceVector(dVector f_vect, iVector row, int num_elem){
	int		cnt;
	double 	*s_vect;

	if(!(s_vect = (dVector) arr_alloc(num_elem+1, double)))
		errMsg("allocation", "reduceVector", "s_vect", 1);

	for (cnt = 1; cnt <= num_elem; cnt++)
		s_vect[cnt] = f_vect[row[cnt]];
	s_vect[0] = oneNorm(s_vect+1, num_elem);

	return s_vect;
}//END reduceVector()

dVector expandVector(dVector red, iVector col, int redElems, int expElems){
	int 	n;
	dVector 	exp;

	if (!(exp = (dVector) arr_alloc(expElems+1, double)) )
		errMsg("allocation", "expandVector", "expanded dVector", 0);

	for (n = 1; n <= redElems; n++ )
		exp[col[n]] = red[n];
	exp[0] = oneNorm(exp+1, expElems);

	return exp;
}//END expandVector

bool equalVector(dVector a, dVector b, int len, double tolerance) {
	int		cnt;

	for (cnt = 1; cnt <= len; cnt++)
		if ( DBL_ABS(a[cnt] - b[cnt]) > tolerance )
			return false;

	return true;
}//END equalVector()

bool equalIntvec(iVector a, iVector b, int len) {
	int		cnt;

	for (cnt = 1; cnt <= len; cnt++)
		if ( a[cnt] != b[cnt] )
			return false;

	return true;
}//END equalIntvec()

bool equalLongIntvec(unsigned long *a, unsigned long *b, int len) {
	int		cnt;

	for (cnt = 1; cnt <= len; cnt++)
		if ( a[cnt] != b[cnt] )
			return false;

	return true;
}//END equalLongIntvec()

bool isZeroVector(dVector a, int len, double tolerance) {
	int		cnt;

	for (cnt = 0; cnt < len; cnt++) {
		if ( DBL_ABS(a[cnt]) >= tolerance )
			return false;
		//		else
		//			a[cnt] = 0.0;
	}

	return true;
}//END isZeroVector()

/*This function will check if a dVector is integer with a predefined gap */
bool isInteger(dVector x, int length, int startIdx, int endIdx, double tolerance){
	int i;

	for (i = startIdx+1; i < endIdx; i++)
		if (fabs(x[i] - round(x[i])) > tolerance)
			return false;

	return true;
}//END isInteger()


dVector duplicVector(dVector a, int len) {
	int		i;
	dVector	b;

	if ((b = (dVector) arr_alloc(len, double))) {
		for (i = 0; i < len; i++)
			b[i] = a[i];
	}

	return b;
}//END duplicArray()

iVector duplicIntvec(iVector a, int len) {
	int		i;
	iVector	b;

	if ((b = (iVector) arr_alloc(len, int))) {
		for (i = 0; i < len; i++)
			b[i] = a[i];
	}

	return b;
}//END duplicArray()

void copyVector(dVector a, dVector b, int len){
	int n;

	for ( n = 0; n < len; n++ )
		b[n] = a[n];

}//END copyVector()

void copyIntvec (iVector a, iVector b, int len) {
	int n;

	for ( n = 0; n < len; n++ )
		b[n] = a[n];

}//END copyVector()

void addVectors(dVector a, dVector b, iVector indices, int len) {
	int n;

	if ( indices == NULL ) {
		for ( n = 1; n <= len; n++ )
			a[n] += b[n];
		a[0] = oneNorm(a+1, len);
	}
	else {
		for ( n = 1; n <= len; n++ )
			a[indices[n]] += b[n];
	}

}//END copy_arr()

void printVector(dVector vec, int len, FILE *fptr){
	int n;

	if ( fptr == NULL ) {
		for ( n = 1; n <= len; n++ )
			printf("%4.6lf, ", vec[n]);
		printf("\n");
	}
	else {
		for ( n = 1; n <= len; n++ )
			fprintf(fptr, "%4.6lf, ", vec[n]);
		fprintf(fptr, "\n");
	}

}//END printVector()

void printVectorWName(dVector vec, cString *vecName, int len, FILE *fptr){
	int n;

	for ( n = 1; n <= len; n++ ) {
		fprintf(fptr, "%s\t\t%4.6lf\n ", vecName[n-1],vec[n]);
		fprintf(fptr, "\n");
	}

}//END printVectorWName()

void printIntvec(iVector vec, int len, FILE *fptr){
	int n;

	if (fptr == NULL) {
		for ( n = 1; n <= len; n++ )
			printf("%d, ", vec[n]);
		printf("\n");
	}
	else {
		for ( n = 1; n <= len; n++ )
			fprintf(fptr, "%d, ", vec[n]);
		fprintf(fptr, "\n");
	}

}//END printIntvec()

void printSparseVector(dVector vec, iVector indices, int len) {
	int n;

	for ( n = 1; n <= len; n++ )
		printf("%4.3lf\t", vec[indices[n]]);
	printf("\n");

}//END printSparseVector()

void printSparseMatrix(sparseMatrix *V, char *cString) {
	int 	cnt;
	printf("%s (%d) ::\n\t\n", cString, V->cnt);
	for (cnt = 1; cnt <= V->cnt; cnt++){
		printf("(%d, %d, %.2f)\t\n", V->row[cnt], V->col[cnt], V->val[cnt]);
		if ( cnt % 5 == 0 )
			printf("\n");
	}
	printf("\n");
}// END print_sparseMatrix()

void printLine() {

	printf("-------------------------------------------------------------------------- \n");

}//END printLine

/* Create a vector on non-zero elements in the input integer vector */
iVector findElements(iVector allElem, int totalElem, int *numUniq){
	iVector	elemUniq;
	int		n, m, len;

	if(!(elemUniq = arr_alloc(totalElem+1, int)))
		errMsg("allocation", "find_cols", "colUniq", 1);

	len = 0;
	/* Copy over all the distinct non-zero elements of allElem */
	for ( n = 1; n <= totalElem; n++ ) {
		if ( allElem[n] > 0 ) {
			m = 1;
			while ( m <= len ) {
				if ( allElem[n] == elemUniq[m] )
					break;
				m++;
			}
			if ( m == len+1 )
				elemUniq[++len] = allElem[n];
		}
	}

	/* Shrink the array down to the number of distinct elements found */
	elemUniq = (iVector) mem_realloc(elemUniq, (len+1)*sizeof(int));

	elemUniq[0] = 0;
	*numUniq = len;

	return elemUniq;
}//END findElems()

/* check if an integer is in an element of input iVector */
bool isInVec(iVector vec, int len, int val) {
	int i;

	for (i = 0; i < len; i++)
		if (vec[i] == val) {
			return true;
		}

	return false;
}//END isInVec()

/* The function encodes an integer dVector _stream_ of given length _len_ into an unsigned long dVector _codeWord_. The _maxValue indicates the maximum value of decoded integer.*/
unsigned long *encodeIntvec(iVector stream, int len, int wordLength, int maxValue) {
	unsigned long *codeWord, temp;
	int j, group, shift, codeLength, numBits;

	numBits = (int) ceil(log2(maxValue));
	codeLength = ceil((double) numBits*len/ (double) wordLength) + 1;

	if ( !(codeWord = (unsigned long *) arr_alloc(codeLength, unsigned long)))
		errMsg("allocation", "encodeIntVec", "codeWord", 0);

	for (j = 1; j <= len; j++) {
		group = numBits*(j-1)/wordLength + 1;
		shift = wordLength - numBits*(j % wordLength);
		temp = (unsigned long) stream[j] << shift;
		codeWord[group] |= temp;
	}

	return codeWord;
}//END encodeIntvec()

/* Decode the column and return the total number of 1's */
iVector decodeIntvec(unsigned long *codeWord, int len, int wordLength, int maxValue) {
	iVector 	stream;
	int 	j, group, shift, numBits, mask = 0;

	numBits = (int) ceil(log2(maxValue));
	for ( j = 0; j < numBits; j++ )
		mask = (mask << 1) + 1;
	stream = (iVector) arr_alloc(len+1, int);

	/* Let's decode phi_col */
	for (j = 1; j <= len; j++) {
		group = numBits*(j-1)/wordLength + 1;
		shift = wordLength - numBits*(j % wordLength);
		stream[j] = (unsigned long) codeWord[group] >> shift;
		stream[j] = stream[j] & mask;
	}

	return stream;
}//END decodeIntvec()

/* This subroutine extracts elements which are common to the two input integer dVectors _a_ and _b_ */
iVector iVectorIntersect(iVector a, iVector b, int lenA, int lenB) {
	iVector inter;
	int	cnt, n;

	if ( !(inter = (iVector) arr_alloc(maximum(lenA, lenB)+1, int)) )
		errMsg("allocation", "iVectorIntersect", "inter", 0);

	cnt = 1;
	for ( n = 1; n <= lenA; n++ )
		if ( isElementIntvec(b, lenB, a[n]) )
			inter[cnt++] = a[n];

	return inter;

}//END iVectorIntersect()

/* This subroutine checks to see if a integer scalar is an element of integer dVector. If so, the subroutine will return the index. If not, a value of -1 is returned. */
int isElementIntvec(iVector vec, int lenVec, int elem) {
	int n = 1;

	while ( n <= lenVec ) {
		if ( vec[n] == elem )
			break;
		n++;
	}

	if ( n > lenVec )
		return -1;
	else
		return n;

}//END isElementIntvec()

void subVectors(dVector a, dVector b, iVector indices, int len){
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

}//END subVectors()


void freeSparseMatrix(sparseMatrix *M) {

	if (M->col) mem_free(M->col);
	if (M->row) mem_free(M->row);
	if (M->val) mem_free(M->val);
	mem_free(M);

}//END freeSparseMatrix()

void freeSparseVector(sparseVector *v) {

	if (v->col) mem_free(v->col);
	if (v->val) mem_free(v->val);
	mem_free(v);

}//END freeSparseMatrix()
