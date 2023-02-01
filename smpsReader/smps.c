/*
 * input.c
 *
 *  Created on: Sep 29, 2015
 *      Author: gjharsha
 */

#include <utils.h>
#include <smps.h>
#include <solver_cplex.h>

int readFiles(cString inputDir, cString probName, oneProblem **orig, timeType **tim, stocType **stoc) {

	/* read problem core file */
	(*orig) = readCore(inputDir, probName);
	if ( (*orig) == NULL ) {
		errMsg("read", "readFiles", "failed to read problem core file", 0);
		return 1;
	}

	/* read problem time file */
	(*tim) = readTime(inputDir, probName, (*orig));
	if ( (*tim) == NULL ) {
		errMsg("read", "readFiles", "failed to read problem time file", 0);
		return 1;
	}

	(*stoc) = readStoc(inputDir, probName, (*orig), (*tim));

#ifdef INPUT_CHECK
	writeStocType((*stoc));
#endif

	return 0;
}//END readFiles()

oneProblem *readCore(cString inputDir, cString probName) {
	LPptr 			lp = NULL;
	char 			probpath[BLOCKSIZE], line[BLOCKSIZE], field1[NAMESIZE], field2[NAMESIZE];
	oneProblem      *orig;
	int				c, nzcnt;
	FILE			*fptr;

	/* Locate the problem core file */
	sprintf(probpath, "%s%s/%s.lp", inputDir, probName, probName);
	fptr = fopen(probpath, "r");
	if ( fptr == NULL ) {
		sprintf(probpath, "%s%s/%s.cor", inputDir, probName, probName);
		fptr = fopen(probpath, "r");
		if ( fptr == NULL ) {
			sprintf(probpath, "%s%s/%s.mps", inputDir, probName, probName);
			fptr = fopen(probpath, "r");
			if ( fptr == NULL ) {
				errMsg("read", "readCore", "failed to open problem core file", 0);
				return NULL;
			}
		}
	}

	/* NAME section: read problem name and compare with that read earlier */
	if ( fgets(line, sizeof line, fptr) != NULL )
		sscanf(line, "%s %s", field1, field2);
	else {
		errMsg("read", "readCore", "failed to read the problem name", 0);
		return NULL;
	}
	if ( !(strncmp(field1, "NAME", 4)) )
		if( strcmp(probName, field2) ) {
			errMsg("read", "readCore", "problem name does not match, in NAME section", 0);
			return NULL;
		}
	fclose (fptr);

	openSolver();
	/* Create LP pointer */
	if ((createProblem(probName, &lp))) {
		errMsg("solver", "readCore", "failed to create problem in solver.\n", 0);
		return NULL;
	}

	if ((readProblem(probpath, lp))) {
		errMsg("solver", "readCore", "failed to read and copy the problem data", 0);
		return NULL;
	}

	/* Allocate memory to the elements of problem and assign default values*/
	orig = (oneProblem *) mem_malloc(sizeof(oneProblem));
	orig->lp = lp;
	orig->name = (cString) mem_calloc(NAMESIZE, sizeof(char));

	/* obtain type of problem read */
	orig->type = getProbType(lp);

	/* obtain the problem elements */
	/* (1) problem name */
	strcpy(orig->name, probName);

	/* (2) objective sense */
	orig->objsen = getObjSen(lp);
	if (!(orig->objsen))
		errMsg("solver", "readCore", "failed to obtain the objective sense", 1);

	/* (3) number of rows */
	orig->mar = getNumRows(lp);
	if (!(orig->mar))
		errMsg("solver", "readCore", "failed to obtain the number of rows in the problem", 1);

	/* (4) number of columns */
	orig->mac = getNumCols(lp);
	if (!(orig->mac))
		errMsg("solver", "readCore", "failed to obtain the number of columns in the problem", 1);

	/* (5) number of non-zeros */
	nzcnt = getNumnz(lp);

	/* continue allocating memory to the elements of problem and assign default values*/
	orig->objx = (dVector) mem_calloc(orig->mac, sizeof(double));
	orig->rhsx = (dVector) mem_calloc(orig->mar, sizeof(double));
	orig->senx = (cString) mem_malloc(orig->mar*sizeof(char));
	orig->matbeg = (iVector) mem_malloc(orig->mac*sizeof(int));
	orig->matcnt = (iVector) mem_malloc(orig->mac*sizeof(int));
	orig->matind = (iVector) mem_malloc(nzcnt*sizeof(int));
	orig->matval = (dVector) mem_malloc(nzcnt*sizeof(double));
	orig->bdl = (dVector) mem_malloc(orig->mac*sizeof(double));
	orig->bdu = (dVector) mem_malloc(orig->mac*sizeof(double));
	orig->ctype = (cString) mem_malloc(orig->mac*sizeof(char));

	/* (6) objective function coefficients */
	if ( (getObjx(lp, 0, orig->mac, orig->objx)) )
		errMsg("solver", "readCore", "failed to obtain the objective coefficients", 1);

	/* (7) constraint right hand side */
	if ( (getRhsx(lp, 0, orig->mar, orig->rhsx)) )
		errMsg("solver", "readCore", "failed to obtain problem constraint right hand sides", 1);

	/* (8) constraint sense */
	if ( (getSense(lp, 0, orig->mar, orig->senx)))
		errMsg("solver", "readCore", "failed to obtain problem constraint sense", 1);

	/* (10) constraint matrix coefficients */
	if ( (getCols(lp, 0, orig->mac, orig->matbeg, orig->matind, orig->matval, nzcnt)) )
		errMsg("solver", "readCore", "failed to obtain the constraint matrix coefficients", 1);
	for ( c = 0; c < orig->mac - 1; c++ )
		orig->matcnt[c] = orig->matbeg[c+1] - orig->matbeg[c];
	orig->matcnt[c] = nzcnt - orig->matbeg[c];

	/* (11) problem variable bounds */
	if( getLb(lp, 0, orig->mac, orig->bdl) )
		errMsg("solver", "readCore", "failed to obtain the problem lower bounds", 1);
	if( getUb(lp, 0, orig->mac, orig->bdu) )
		errMsg("solver", "readCore", "failed to obtain the problem upper bounds", 1);

	if ( orig->type == PROB_MILP || orig->type == PROB_MIQP ) {
		/* (12) get problem constraint type */
		if ( getCtype(lp, 0, orig->mac, orig->ctype) )
			errMsg("solver", "readCore", "failed to obtain the variable types", 1);
		orig->numInt = getNumInt(lp);
		orig->numBin = getNumBinary(lp);
	}
	else {
		for ( c = 0; c < orig->mac; c++ )
			orig->ctype[c] = 'C';
		orig->numInt = 0;
	}

	/* Allocate memory to hold the names of problem elements */
	orig->objname = (cString) mem_calloc(NAMESIZE, sizeof(char));
	orig->cstorsz = -getCstoreSize(lp, 0, orig->mac);
	if ( orig->cstorsz <= 0 )
		errMsg("solver", "readCore", "Could not determine amount of space for column names", 1);
	orig->cname = (cString *) mem_malloc(orig->mac*sizeof(char *));
	orig->cstore = (cString) mem_malloc(orig->cstorsz);

	orig->rstorsz = -getRstoreSize(lp, 0, orig->mar);
	if ( orig->rstorsz < 0 )
		errMsg("solver", "readCore", "Could not determine amount of space for row names", 1);
	orig->rname = (cString *) mem_malloc(orig->mar*sizeof(char *));
	orig->rstore = (cString) mem_malloc(orig->rstorsz);

	/* (12) objective name */
	if ( (getObjName(lp, orig->objname)) )
		errMsg("solver", "readCore", "failed to obtain objective name", 1);

	/* (13) problem row name */
	if ( (getRowName(lp, 0, orig->mar, orig->rname, orig->rstore, orig->rstorsz)) )
		errMsg( "solver", "readCore", "failed to obtain row names", 1);

	/* (14) problem column name */
	if ( (getColName(lp, 0, orig->mac, orig->cname, orig->cstore, orig->cstorsz)) )
		errMsg("solver", "readCore", "failed to obtain column names", 1);

	orig->matsz = nzcnt;
	orig->macsz = nzcnt;
	orig->marsz = nzcnt;
	orig->numnz = nzcnt;

	return orig;
}//END readCore()

timeType *readTime(cString inputDir, cString probName, oneProblem *orig) {
	timeType	*tim;
	char		probpath[2*BLOCKSIZE], line[BLOCKSIZE], field1[NAMESIZE], field2[NAMESIZE];
	int			defaultStages = 20, n, m;
	FILE		*fptr;

	/* Locate the problem core file */
	sprintf(probpath, "%s%s/%s.tim", inputDir, probName, probName);

	/* open the time file */
	fptr = fopen(probpath,"r");
	if (fptr == NULL) {
		errMsg("read","readTime","failed to read time file/missing file", 0);
		return NULL;
	}

	/* allocate memory and initialize */
	if (!(tim = (timeType *) mem_malloc(sizeof(timeType))))
		errMsg("allocation", "readTime", "timeType",0);
	if(!(tim->stgNames = (cString *) mem_malloc(defaultStages*sizeof(cString))))
		errMsg("allocation", "readTime", "stgNames in timeType", 0);
	tim->numStages = 0; n = 0;
	tim->numCols = 0; tim->numRows = 0;

	/* TIME section: collect the problem name */
	if ( fgets(line, sizeof line, fptr) != NULL )
		sscanf(line, "%s %s", field1, field2);
	else {
		errMsg("read", "readTime", "failed to read the problem name", 0);
		return NULL;
	}
	if ( !(strncmp(field1, "TIME", 4)) )
		if( strcmp(probName, field2) ) {
			errMsg("read", "readTime", "problem name does not match, TIME section", 0);
			return NULL;
		}

	/* PERIODS section: collect the time file type */
	if ( fgets(line, sizeof line, fptr) != NULL )
		sscanf(line, "%s %s", field1, field2);
	else {
		errMsg("read", "readTime", "failed to read the time file type", 0);
		return NULL;
	}

	if ( !(strncmp(field1, "PERIODS", 7)) ) {
		if ( !(strncmp(field2, "EXPLICIT", 8)) )
			tim->type = 1;
		else
			tim->type = 0;
	}
	else {
		errMsg("read", "readTime", "unknown header record time file, PERIODS section", 0);
		return NULL;
	}

	if ( tim->type == 0 ){
		if( !(tim->row = (iVector) arr_alloc(defaultStages, int)) )
			errMsg("allocation", "readTime", "rowNames in timeType", 0);
		if( !(tim->col = (iVector) arr_alloc(defaultStages, int)) )
			errMsg("allocation", "readTime", "colNames in timeType", 0);
		while ( fgets(line, sizeof line, fptr )!= NULL ) {
			if (line[0] != '*' && strncmp(line,"ENDATA",6)) /* If it is not a comment line and end of data in the file proceed to read the contents */ {
				if ( !(tim->stgNames[n] =  (cString) mem_malloc(NAMESIZE*sizeof(char))))
					errMsg("allocation", "readTime", "individual stage names", 0);
				sscanf(line, "%s %s %s", field1, field2, tim->stgNames[n]);
				/* find the column and row coordinates in original problem */
				m = 0;
				while ( m < orig->mac ){
					if ( !(strcmp(field1, orig->cname[m])) )
						break;
					m++;
				}
				if ( m == orig->mac ) {
					errMsg("read", "readTime", "unknown column name in the time file", 0);
					return NULL;
				}
				tim->col[n] = m;
				if ( !(strcmp(field2, orig->objname)) )
					m = 0;
				else {
					m = 0;
					while (m < orig->mar ) {
						if ( !(strcmp(field2, orig->rname[m])) )
							break;
						m++;
					}
				}
				if ( m == orig->mar ) {
					errMsg("read", "readTime", "unknown row name in the time file", 0);
					return NULL;
				}
				tim->row[n] = m;

				/* increment stage counter and proceed */
				tim->numStages++; n++;
				if ( n == defaultStages ) {
					errMsg("allocation", "readTime",
							"ran out of memory to store periods in implicit declaration; increase default value", 0);
					return NULL;
				}
			}
		}
		tim->numCols = tim->numRows = tim->numStages;
		tim->rowStg = NULL; tim->colStg = NULL;
	}
	else {
		errMsg("read", "readTime", "explicit time file description is currently not supported", 0);
		/* TODO (HG): complete the code to read explicit time file declaration */
		return NULL;
	}

	/* reallocate memory elements of time structure */
	tim->stgNames = (cString *) mem_realloc(tim->stgNames, tim->numStages*sizeof(cString));
	if (tim->type == 0) {
		tim->col = (iVector) mem_realloc(tim->col, tim->numStages*sizeof(int));
		tim->row = (iVector) mem_realloc(tim->row, tim->numStages*sizeof(int));
	}
	else {
		tim->col = (iVector) mem_realloc(tim->col, tim->numStages*sizeof(int));
		tim->row = (iVector) mem_realloc(tim->row, tim->numStages*sizeof(int));
		tim->colStg = (iVector) mem_realloc(tim->colStg, tim->numCols*sizeof(int));
		tim->rowStg = (iVector) mem_realloc(tim->rowStg, tim->numRows*sizeof(int));
	}
	fclose(fptr);
	return tim;
}//END readTime()

stocType *readStoc(cString inputDir, cString probName, oneProblem *orig, timeType *tim) {
	stocType *stoc;
	cString 	*rvRows = NULL, *rvCols = NULL, *fields = NULL;
	char	probpath[2*BLOCKSIZE], line[BLOCKSIZE], fieldType;
	FILE	*fptr;
	int		maxOmegas = 1000, maxVals = 4000, n, numFields, maxFields = 10;

	/* Locate the problem sto file */
	sprintf(probpath, "%s%s/%s.sto", inputDir, probName, probName);

	/* open the time file */
	fptr = fopen(probpath,"r");
	if (fptr == NULL) {
		errMsg("read","readStoch","failed to read stoch file file/missing file", 0);
		return NULL;
	}

	/* allocate memory to field locations */
	if (!(fields = (cString *) arr_alloc(maxFields, cString)) )
		errMsg("allocation", "readStoc", "field locations", 0);
	for (n = 0; n < maxFields; n++ )
		if ( !(fields[n] = (cString) arr_alloc(NAMESIZE, char)) )
			errMsg("allocation", "readStoc", "individual field location", 0);

	/* allocate memory to stocType and initialize elements */
	if ( !(stoc = (stocType *) mem_malloc(sizeof(stocType))) )
		errMsg("allocation", "readStoc", "stoc", 0);
	if ( !(stoc->type = (cString) arr_alloc(NAMESIZE, char)) )
		errMsg("allocation", "readStoc", "stoc->type", 0);
	if ( !(stoc->col = (iVector) arr_alloc(maxOmegas, int)) )
		errMsg("allocation", "readStoc", "stoc->col", 0);
	if ( !(stoc->row = (iVector) arr_alloc(maxOmegas, int)) )
		errMsg("allocation", "readStoc", "stoc->row", 0);
	if ( !(stoc->mean = (dVector) arr_alloc(maxOmegas, double)) )
		errMsg("allocation", "readStoc", "stoc->mean", 0);
	if ( !(stoc->groupBeg = (iVector) arr_alloc(maxFields, int)) )
		errMsg("allocation", "readStoc", "stoc->groupBeg", 0);
	if ( !(stoc->numPerGroup = (iVector) arr_alloc(maxFields, int)) )
		errMsg("allocation", "readStoc", "stoc->numPerGroup", 0);
	stoc->numOmega = 0;
	stoc->numGroups = 0;
	stoc->sim = false;
	stoc->isDiscrete = true;
	stoc->mod = NULL;

	/* STOCH section: read problem name and compare with that read earlier */
	if ( fgets(line, sizeof line, fptr) != NULL )
		sscanf(line, "%s %s", fields[0], fields[1]);
	else {
		errMsg("read", "readStoc", "failed to read the problem name", 0);
		return NULL;
	}
	if ( !(strncmp(fields[0], "STOCH", 5)) )
		if( strcmp(probName, fields[1]) ) {
			errMsg("read", "readStoc", "Problem name do not match, in STOCH section", 0);
			return NULL;
		}

	while ( !(getLine(&fptr, fields, &fieldType, &numFields)) ) {
		START_OVER: /* Used to continue parsing the stoch file when there are many stochastic elements. */
		if ( !(strcmp(fields[0], "INDEP")) ) {
			if ( readIndep(fptr, fields, orig, maxOmegas, maxVals, stoc, &rvRows, &rvCols) ) {
				errMsg("read", "readStoc", "failed to read stoch file with independent data", 0);
				return NULL;
			}
			goto START_OVER;
		}
		if ( !(strcmp(fields[0], "BLOCKS")) ) {
			if ( readBlocks(fptr, fields, orig, maxOmegas, maxVals, stoc, &rvRows, &rvCols) ) {
				errMsg("read", "readStoc", "failed to read stoch file with blocks", 0);
				return NULL;
			}
			goto START_OVER;
		}
		if ( !(strcmp(fields[0], "SCENARIOS")) ) {
			if ( readScenarios(fptr, fields, orig, tim, maxOmegas, maxVals, stoc) ) {
				errMsg("read", "readStoc", "failed to read stoch file with scenarios", 0);
				return NULL;
			}
			goto START_OVER;
		}
		if ( !(strcmp(fields[0], "ENDATA")) )
			break;
	}

	/* free allocated memory */
	maxOmegas = stoc->numOmega;
	if ( stoc->mod != NULL )
		maxOmegas += stoc->mod->M;
	if ( rvCols ) {
		for ( n = 0; n < maxOmegas; n++ ) {
			if(rvCols[n]) mem_free(rvCols[n]);
		}
		mem_free(rvCols);
	}
	if ( rvRows ) {
		for ( n = 0; n < maxOmegas; n++ ) {
			if(rvRows[n]) mem_free(rvRows[n]);
		}
		mem_free(rvRows);
	}
	if ( fields ) {
		for (n = 0; n < maxFields; n++ )
			if (fields[n]) mem_free(fields[n]);
		mem_free(fields);
	}
	fclose(fptr);

	/* Reallocate memory to elements of stocType */
	stoc->row  = (iVector) mem_realloc(stoc->row, stoc->numOmega*sizeof(int));
	stoc->col  = (iVector) mem_realloc(stoc->col, stoc->numOmega*sizeof(int));
	stoc->mean = (dVector) mem_realloc(stoc->mean, stoc->numOmega*sizeof(double));

	stoc->numPerGroup = (iVector) mem_realloc(stoc->numPerGroup, stoc->numGroups*sizeof(int));
	stoc->groupBeg 	  = (iVector) mem_realloc(stoc->groupBeg ,stoc->numGroups*sizeof(int));

	return stoc;
}//END readStoc()

int readIndep(FILE *fptr, cString *fields, oneProblem *orig, int maxOmegas, int maxVals, stocType *stoc, cString 	**rvRows, cString **rvCols) {

	/* Mark where the group beings */
	stoc->groupBeg[stoc->numGroups] = stoc->numOmega;

	/* allocate memory to hold the names of random variable */
	if ( !((*rvRows) = (cString *) arr_alloc(maxOmegas, cString)) )
		errMsg("allocation", "readIndep", "rvNames", 0);
	if ( !((*rvCols) = (cString *) arr_alloc(maxOmegas, cString)) )
		errMsg("allocation", "readIndep", "rvNames", 0);

	if ( !(strcmp(fields[1], "DISCRETE")) ) {
		if ( readIndepDiscrete(fptr, fields, maxOmegas, maxVals, rvRows, rvCols, orig, stoc)) {
			errMsg("read", "readIndep", "failed to read independent discrete random variables", 0);
			return 1;
		}
	}
	else if ( strstr(fields[1], "NORMAL") != NULL ) {
		if( readNormal(fptr, fields, maxOmegas, rvRows, rvCols, orig, stoc) ) {
			errMsg("read", "readIndep", "failed to read independent normal distribution", 0);
			return 1;
		}
	}
	else if ( !(strcmp(fields[1], "EXPONENTIAL")) ) {
		errMsg("read", "readIndeps", "no support for exponential distribution type in INDEP section", 0);
		return 1;
	}
	else if ( !(strcmp(fields[1], "UNIFORM")) ) {
		errMsg("read", "readIndeps", "no support for uniform distribution type in INDEP section", 1);
		return 1;
	}
	else if ( !(strcmp(fields[1], "GAMMA")) ) {
		errMsg("read", "readIndeps", "no support for gamma distribution type in INDEP section", 1);
		return 1;
	}
	else if ( !(strcmp(fields[1], "GEOMETRIC")) ) {
		errMsg("read", "readIndeps", "no support for geometric distribution type in INDEP section", 1);
		return 1;
	}
	else {
		errMsg("read", "readIndeps", "unknown distribution type in INDEP section", 1);
		return 1;
	}

	/* increase the number of stochastic variables groups */
	stoc->numPerGroup[stoc->numGroups] = stoc->numOmega - stoc->groupBeg[stoc->numGroups];
	stoc->numGroups++;

	return 0;
}//END readIndep()

int readIndepDiscrete(FILE *fptr, cString *fields, int maxOmegas, int maxVals, cString **rvRows, cString **rvCols, oneProblem *orig, stocType *stoc) {
	int numFields, n;
	char strType;

	/* store the type of stochastic process encountered */
	sprintf(stoc->type, "INDEP_DISCRETE");
	stoc->sim = false;

	stoc->numVals = (iVector) arr_alloc(maxOmegas, int);
	stoc->vals    = (dVector *) arr_alloc(maxOmegas, dVector);
	stoc->probs   = (dVector *) arr_alloc(maxOmegas, dVector);
	stoc->mod = NULL;

	while (true) {
		getLine(&fptr, fields, &strType, &numFields);
		if (strType != 'f')
			break;										// Encountered ENDATA or a new group of random variables
		n = stoc->numOmega - 1;
		if ( n > maxOmegas ) {
			errMsg("allocation", "readIndep", "reached maxOmega limit for INDEP format", 0);
			return 1;
		}
		while (n >= 0 ) {
			if ( !(strcmp(fields[0], (*rvCols)[n])) && !(strcmp(fields[1], (*rvRows)[n])) )
				break;
			n--;
		}
		if ( n == -1 ) {
			/* new random variable encountered */
			if ( !((*rvRows)[stoc->numOmega] = (cString) arr_alloc(NAMESIZE, char)) )
				errMsg("allocation", "readIndep", "rvNames[n]", 0);
			if ( !((*rvCols)[stoc->numOmega] = (cString) arr_alloc(NAMESIZE, char)) )
				errMsg("allocation", "readIndep", "rvNames[n]", 0);
			if ( !(stoc->vals[stoc->numOmega] = (dVector) arr_alloc(maxVals, double)) )
				errMsg("allocation", "readIndep","omega.vals[n]", 0);
			if ( !(stoc->probs[stoc->numOmega] = (dVector) arr_alloc(maxVals, double)) )
				errMsg("allocation", "readIndep", "omega.probs[n]", 0);

			strcpy((*rvCols)[stoc->numOmega], fields[0]);
			strcpy((*rvRows)[stoc->numOmega], fields[1]);
			stoc->numVals[stoc->numOmega++] = 0;
			/* identify row and column coordinates in the problem */
			if ( !(strcmp(fields[0], "RHS")) )
				n = -1;
			else {
				n = 0;
				while ( n < orig->mac ){
					if ( !(strcmp((*rvCols)[stoc->numOmega-1], orig->cname[n])) )
						break;
					n++;
				}
			}
			if ( n == orig->mac ) {
				errMsg("read", "readIndep", "unknown column name in the stoch file", 0);
				return 1;
			}
			stoc->col[stoc->numOmega-1] = n;
			if ( !(strcmp(fields[1], orig->objname)) )
				n = -1;
			else {
				n = 0;
				while (n < orig->mar ) {
					if ( !(strcmp((*rvRows)[stoc->numOmega-1], orig->rname[n])) )
						break;
					n++;
				}
			}
			if ( n == orig->mar ) {
				errMsg("read", "readIndep", "unknown row name in the stoch file", 0);
				return 1;
			}
			stoc->row[stoc->numOmega-1] = n;
		}
		if ( numFields == 4) {
			stoc->vals[stoc->numOmega-1][stoc->numVals[stoc->numOmega-1]] = str2float(fields[2]);
			stoc->probs[stoc->numOmega-1][stoc->numVals[stoc->numOmega-1]] = str2float(fields[3]);
			stoc->mean[stoc->numOmega-1] += str2float(fields[2])*str2float(fields[3]);
			stoc->numVals[stoc->numOmega-1]++;
		}
		else if ( numFields == 5 ) {
			stoc->vals[stoc->numOmega-1][stoc->numVals[stoc->numOmega-1]] = str2float(fields[2]);
			stoc->probs[stoc->numOmega-1][stoc->numVals[stoc->numOmega-1]] = str2float(fields[4]);
			stoc->mean[stoc->numOmega-1] += str2float(fields[2])*str2float(fields[4]);
			stoc->numVals[stoc->numOmega-1]++;
		}
		else {
			errMsg("read", "readIndep", "missing field in stoch file", 0);
			return 1;
		}
	}

	/* Reallocate memory to fit the exact size */
	stoc->numVals = (iVector) mem_realloc(stoc->numVals, stoc->numOmega*sizeof(int));
	for ( n = 0; n < stoc->numOmega; n++ ) {
		stoc->vals[n]  = (dVector) mem_realloc(stoc->vals[n], stoc->numVals[n]*sizeof(double));
		stoc->probs[n] = (dVector) mem_realloc(stoc->probs[n], stoc->numVals[n]*sizeof(double));
	}
	stoc->vals  = (dVector *) mem_realloc(stoc->vals, stoc->numOmega*sizeof(dVector));
	stoc->probs = (dVector *) mem_realloc(stoc->probs, stoc->numOmega*sizeof(dVector));

	return 0;
}//END readIndepDiscrete()

int readNormal(FILE *fptr, cString *fields, int maxOmegas, cString **rvRows, cString **rvCols, oneProblem *orig, stocType *stoc) {
	int n, numFields;
	char strType;

	/* continuous distribution, use a simulator */
	stoc->sim = false;
	stoc->isDiscrete = false;
	sprintf(stoc->type, "INDEP_%s",fields[1]);

	stoc->vals	  = (dVector *) arr_alloc(1, dVector);
	stoc->vals[0] = (dVector) arr_alloc(maxOmegas, double);

	stoc->probs = NULL; stoc->numVals = NULL; stoc->mod = NULL;

	while (true) {
		getLine(&fptr, fields, &strType, &numFields);
		if (strType != 'f')
			break;
		n = stoc->numOmega - 1;
		if ( n > maxOmegas ) {
			errMsg("allocation", "readIndep", "reached maxOmega limit for INDEP format", 0);
			return 1;
		}
		while (n >= 0 ) {
			if ( !(strcmp(fields[0], (*rvCols)[n])) && !(strcmp(fields[1], (*rvRows)[n])) )
				break;
			n--;
		}
		if ( n == -1 ) {
			/* new random variable encountered */
			if ( stoc->numOmega == maxOmegas ) {
				errMsg("read", "readIndep", "ran out of memory to store row and column names", 0);
				return 1;
			}
			if ( !((*rvRows)[stoc->numOmega] = (cString) arr_alloc(NAMESIZE, char)) )
				errMsg("allocation", "readIndep", "rvNames[n]", 0);
			if ( !((*rvCols)[stoc->numOmega] = (cString) arr_alloc(NAMESIZE, char)) )
				errMsg("allocation", "readIndep", "rvNames[n]", 0);

			strcpy((*rvCols)[stoc->numOmega], fields[0]);
			strcpy((*rvRows)[stoc->numOmega], fields[1]);
			stoc->numOmega++;

			/* Check to see if the random variable corresponds to error terms in a linear transformation or ARMA model. */
			if ( strcmp(fields[1], "LAGGED") ) {
				/* Identify row and column coordinates in the problem */
				if ( !(strcmp(fields[0], "RHS")) )
					n = -1;
				else {
					n = 0;
					while ( n < orig->mac ){
						if ( !(strcmp((*rvCols)[stoc->numOmega-1], orig->cname[n])) )
							break;
						n++;
					}
				}
				if ( n == orig->mac ) {
					errMsg("read", "readIndep", "unknown column name in the stoch file", 0);
					return 1;
				}
				stoc->col[stoc->numOmega-1] = n;
				if ( !(strcmp(fields[1], orig->objname)) )
					n = -1;
				else {
					n = 0;
					while (n < orig->mar ) {
						if ( !(strcmp((*rvRows)[stoc->numOmega-1], orig->rname[n])) )
							break;
						n++;
					}
				}
				if ( n == orig->mar ) {
					errMsg("read", "readIndep", "unknown row name in the stoch file", 0);
					return 1;
				}
				stoc->row[stoc->numOmega-1] = n;
			}
			else {
				stoc->row[stoc->numOmega-1] = stoc->col[stoc->numOmega-1] = -1;
			}
		}
		if ( numFields == 4) {
			/* note, standard deviation is held in the first vals field */
			stoc->mean[stoc->numOmega-1] 	= str2float(fields[2]);
			stoc->vals[0][stoc->numOmega-1] = sqrt(str2float(fields[3]));
		}
		else if ( numFields == 5 ) {
			stoc->mean[stoc->numOmega-1] 	= str2float(fields[2]);
			stoc->vals[0][stoc->numOmega-1] = sqrt(str2float(fields[4]));
		}
		else {
			errMsg("read", "readIndep", "missing field in stoch file", 0);
			return 1;
		}
	}

	stoc->vals[0] = (dVector) mem_realloc(stoc->vals[0], stoc->numOmega*sizeof(double));

	return 0;
}//END readNormal()

int readBlocks(FILE *fptr, cString *fields, oneProblem *orig, int maxOmegas, int maxVals, stocType *stoc, cString **rvRows, cString **rvCols) {

	/* Mark where the group beings */
	stoc->groupBeg[stoc->numGroups] = stoc->numOmega;

	if ( !(strcmp(fields[1], "DISCRETE")) ) {
		/* store the type of stochastic process encountered */
		sprintf(stoc->type, "BLOCKS_DISCRETE");
		if ( readOneBlock(fptr, fields, orig, maxOmegas, maxVals, true, stoc) ) {
			errMsg("read", "readBlocks", "failed to read independent blocks structure", 0);
			return 1;
		}
	}
	else if ( !(strcmp(fields[1], "MVNORMAL")) ) {
		errMsg("read", "readBlocks", "no support for multivariate normal distribution type in BLOCKS section", 0);
		return 1;
	}
	else if ( !(strcmp(fields[1], "LINTR")) ) {
		if ( readLinTrans(fptr, fields, orig, stoc, maxOmegas, rvRows, rvCols) ) {
			errMsg("read", "readLinTran", "failed to read linear transformation structure.", 0);
			return 1;
		}

		/* increase the number of stochastic variables groups */
		stoc->numPerGroup[stoc->numGroups] = stoc->numOmega - stoc->groupBeg[stoc->numGroups];
		stoc->numGroups++;
	}
	else {
		errMsg("read", "readBlocks", "unknown distribution type in BLOCKS section", 1);
		return 1;
	}

	return 0;
}//END readBlocks()

int readOneBlock(FILE *fptr, cString *fields, oneProblem *orig, int maxOmegas, int maxVals, bool origRV, stocType *stoc) {
	cString 	*rvRows, *rvCols;
	char 	strType, currBlock[NAMESIZE] = "\0";
	int		numFields, numRV=0, n;
	bool	newBlk;

	/* allocate memory to hold the names of random variable */
	if ( !(rvRows = (cString *) arr_alloc(maxOmegas, cString)) )
		errMsg("allocation", "readOneBlock", "rvNames", 0);
	if ( !(rvCols = (cString *) arr_alloc(maxOmegas, cString)) )
		errMsg("allocation", "readOneBlock", "rvNames", 0);

	for ( n = 0; n < maxOmegas; n++) {
		if ( !(rvCols[n] = (cString) arr_alloc(NAMESIZE, char)) )
			errMsg("allocation", "readOneBlock", "rvCols", 0);
		if ( !(rvRows[n] = (cString) arr_alloc(NAMESIZE, char)) )
			errMsg("allocation", "readOneBlock", "rvRows", 0);
	}

	stoc->numVals = (iVector) arr_alloc(maxOmegas, int);
	stoc->vals    = (dVector *) arr_alloc(maxOmegas, dVector);
	stoc->probs   = (dVector *) arr_alloc(maxOmegas, dVector);
	stoc->mod = NULL;

	while (true) {
		getLine(&fptr, fields, &strType, &numFields);
		if (strType != 'f')
			break;
		if ( !(strcmp(fields[0], "BL")) ) {
			/* new realization of the block */
			if ( strcmp(currBlock, fields[1]) ) {
				/* first encounter with the block, prepare to record names of random variables */
				newBlk = true;
				strcpy(currBlock, fields[1]);
				stoc->groupBeg[stoc->numGroups] = stoc->numOmega;
				stoc->numPerGroup[stoc->numGroups] = numRV = 0;
				if ( !(stoc->probs[stoc->numGroups] = (dVector) arr_alloc(maxVals, double)) )
					errMsg("allocation", "readOneBlock", "stoc->prob[n]", 0);
				stoc->probs[stoc->numGroups][stoc->numVals[stoc->numGroups]++] = str2float(fields[3]);
				stoc->numGroups++;
			}
			else {
				newBlk = false;
				if ( stoc->numVals[stoc->numGroups-1] == maxVals )
					errMsg("allocation", "readOneBlock", "exceeded memory limit on maxVals", 1);
				stoc->probs[stoc->numGroups-1][stoc->numVals[stoc->numGroups-1]++] = str2float(fields[3]);
			}
		}
		else {
			/* read block elements */
			if (newBlk) {
				/* record names of random variables on their first realization */
				strcpy(rvCols[numRV], fields[0]);
				if ( origRV ) {
					/* if the random variables are problem random variables then identify their coordinates in original problem */
					strcpy(rvRows[numRV], fields[1]);
				}
				/* column coordinates */
				if ( !(strcmp(fields[0], "RHS")) )
					n = -1;
				else {
					n = 0;
					while ( n < orig->mac ) {
						if ( !(strcmp(fields[0], orig->cname[n])) )
							break;
						n++;
					}
				}
				stoc->col[stoc->numOmega] = n;
				/* row coordinates */
				if ( !(strcmp(fields[1], orig->objname)) )
					n = -1;
				else {
					n = 0;
					while ( n < orig->mar ) {
						if ( !(strcmp(fields[1], orig->rname[n])) )
							break;
						n++;
					}
				}
				stoc->row[stoc->numOmega] = n;

				/* make sure there is memory space available for new realization and store it */
				if (stoc->numOmega == maxOmegas )
					errMsg("allocation", "readOneBlock", "reached max limit maxOmegas", 1);

				if ( !(stoc->vals[stoc->numOmega] = (dVector) arr_alloc(maxVals, double)) )
					errMsg("allocation", "readOneBlock","omega.vals[n]", 0);

				if (origRV == 1) {
					stoc->vals[stoc->numOmega][stoc->numVals[stoc->numGroups-1]-1] = str2float(fields[2]);
					stoc->mean[stoc->numOmega] += stoc->probs[stoc->numGroups-1][stoc->numVals[stoc->numGroups-1]-1]*str2float(fields[2]);
					stoc->numOmega++;
				}
				else {
					errMsg("read", "readOneBlock", "reading auxiliary variables not supported", 0);
					return 1;
				}
				/* increment the number of random variables in the group */
				stoc->numPerGroup[stoc->numGroups-1]++;
				numRV++;
			}
			else {
				/* locate the random variable in the list and record realization */
				n = 0;
				while (n < numRV) {
					if ( origRV == 	false && !(strcmp(rvCols[n], fields[0])) )
						break;
					else if ( origRV == true && !(strcmp(rvCols[n], fields[0])) && !(strcmp(rvRows[n], fields[1])) )
						break;
					n++;
				}
				if ( n == numRV )
					errMsg("read", "readOneBlock", "unknown block random variable name", 1);

				n += stoc->groupBeg[stoc->numGroups-1];
				if ( origRV == 1 ) {
					/* the third field has values */
					stoc->vals[n][stoc->numVals[stoc->numGroups-1]-1] = str2float(fields[2]);
					stoc->mean[n] += stoc->probs[stoc->numGroups-1][stoc->numVals[stoc->numGroups-1]-1]*str2float(fields[2]);
				}
				else {
					/* the second field has values */
					errMsg("read", "readOneBlock", "reading auxiliary variables not supported", 0);
					return 1;
				}
			}
		}
	}

	for ( n = 0; n < maxOmegas; n++) {
		if ( rvRows[n] ) mem_free(rvRows[n]);
		if ( rvCols[n] ) mem_free(rvCols[n]);
	}
	mem_free(rvRows); mem_free(rvCols);

	return 0;
}//END readOneBlock()

/* The subroutine reads the linear transformation matrix information for stochastic processes of the following form:
 *
 * 				y_t = c_t + \sum_{j=1}^m \phi_j*y_{t-j} + \sum_{j=0}^n \theta_j*\epsilon_{t-j}
 *
 * where, \epsilon_t follows a normal distribution with mean (\mu, \sigma) known as the error/residual random variable.
 *
 * The subroutine assumes that the stoc file begins by first describing the residual random variable. This is stored as
 * the first group of random variables in our stocType structure.
 */
int readLinTrans(FILE *fptr, cString *fields, oneProblem *orig, stocType *stoc, int maxOmegas, cString **rvRows, cString **rvCols) {
	statModel *model;
	iVector 	periodBeg;
	char 	strType, currBlock[NAMESIZE] = "\0", currLag[NAMESIZE] = "\0";
	int		numFields, period, numPeriods = 0, maxP = 10, maxQ = 10,
			maxMatcnt = stoc->numOmega, j, col, row, offset, maxPeriods = 50;
	bool	newLag;

	periodBeg = (iVector) arr_alloc(maxPeriods, int);

	/* allocate memory to hold information about the linear transformation stochastic process */
	if ( !(model = (statModel *) mem_malloc(sizeof(statModel))) )
		errMsg("allocation", "readLinTrans", "statModel", 0);
	model->eta = model->sigma = NULL; model->muEps = NULL; model->cvEps = NULL;
	model->AR = model->MA = NULL;
	model->p = model->q = 0; model->N = model->M = 0;

	offset = stoc->numOmega;
	/* In this case, the previous group of random variable is treated as the residual/noise random variables.
	 * These random variables are moved to the statModel. */
	if ( strstr(stoc->type, "INDEP_NORMAL") != NULL ) {
		model->M = stoc->numPerGroup[stoc->numGroups-1];

		/* Setup the mean dVector and the covariance matrix of residual process/noise */
		model->muEps = (dVector) arr_alloc(model->M, double);
		model->cvEps = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix));
		model->cvEps->col = (iVector) arr_alloc(model->M, int);
		model->cvEps->row = (iVector) arr_alloc(model->M, int);
		model->cvEps->val = (dVector) arr_alloc(model->M, double);

		for ( j = 0; j < model->M; j++ ) {
			model->muEps[j] = stoc->mean[j];
			model->cvEps->col[j] = model->cvEps->row[j] = j; /* Only diagonal elements are non-zero under independence assumption */
			model->cvEps->val[j] = stoc->vals[0][j];
		}

		/* Remove the previous group of random variables are the main list in stocType */
		stoc->numGroups--;
		stoc->numOmega -= stoc->numPerGroup[stoc->numGroups];
		stoc->groupBeg[stoc->numGroups] = stoc->numOmega;
		stoc->vals[0] = (dVector) mem_realloc(stoc->vals[0], maxOmegas*sizeof(dVector));
	}
	else {
		errMsg("read", "readLinTrans", "currently only independent normal residual processes are supported", 0);
		return 1;
	}

	/* Update the stocType */
	strcpy(stoc->type, "LINTRAN");
	stoc->sim = false;

	/* Read from the stoc file line-by-line */
	while (true) {
		getLine(&fptr, fields, &strType, &numFields);
		if (strType != 'f')
			break;
		if ( !(strcmp(fields[0], "BL")) ) {
			/* New block of encountered of type 'BL': random variables in a particular time period/stage (i.e., elements of y_t) */
			strcpy(currBlock, fields[0]);
			periodBeg[numPeriods++] = stoc->numOmega;
			if ( numPeriods == maxPeriods ) {
				maxPeriods *= 2;
				periodBeg = (iVector) mem_realloc(periodBeg, maxPeriods*sizeof(int));
			}
			model->N = stoc->numOmega - model->N;
		}
		else if ( !(strcmp(fields[0], "RV")) || !(strcmp(fields[0], "HV")) || !(strcmp(fields[0], "LV")) )  {
			/* New block of encountered of type
			 * 		'RV': following elements will be for matrix \theta_0
			 * 		'HV': following elements will be for matrix \phi_j
			 * 		'LV': following elements will be for matrix \theta_j		*/
			strcpy(currBlock, fields[0]);
			if ( strcmp(currLag, fields[3]) ) {
				strcpy(currLag, fields[3]);
				newLag = true;
			}
			else
				newLag = false;

			/* Find the random variable to which the column of transformation matrix corresponds to. */
			col = 0;
			while ( (strcmp((*rvCols)[col], fields[1])) || (strcmp((*rvRows)[col], fields[2])) )
				col++;

			if ( col >= offset ) { /* Check to make sure that the column is not a noise random variable */
				col -= model->M;
				period = 0;
				while ( period < numPeriods ) {
					if ( col < periodBeg[period] )
						break;
					period++;
				}
				col -= periodBeg[period-1];
			}
		}
		else {
			/* Block data */
			if ( !(strcmp(currBlock, "BL")) ) {
				/* Read the column and row names and identify their coordinates in the original problem. */
				if ( !(strcmp(fields[0], "RHS")) )
					j = -1;
				else {
					j = 0;
					while ( j < orig->mac ) {
						if ( !(strcmp(fields[0], orig->cname[j])) )
							break;
						j++;
					}
				}
				stoc->col[stoc->numOmega] = j;
				if ( !(strcmp(fields[1], orig->objname)) )
					j = -1;
				else {
					j = 0;
					while ( j < orig->mar ) {
						if ( !(strcmp(fields[1], orig->rname[j])) )
							break;
						j++;
					}
				}
				stoc->row[stoc->numOmega] = j;
				if ( !((*rvRows)[model->M+stoc->numOmega] = (cString) arr_alloc(NAMESIZE, char)) )
					errMsg("allocation", "readIndep", "rvNames[n]", 0);
				strcpy((*rvRows)[model->M+stoc->numOmega], fields[1]);
				if ( !((*rvCols)[model->M+stoc->numOmega] = (cString) arr_alloc(NAMESIZE, char)) )
					errMsg("allocation", "readIndep", "rvNames[n]", 0);
				strcpy((*rvCols)[model->M+stoc->numOmega], fields[0]);
				stoc->vals[0][stoc->numOmega] = str2float(fields[2]);
				stoc->numOmega++;

				/* make sure there is memory space available for new realization and store it */
				if (stoc->numOmega == maxOmegas )
					errMsg("allocation", "readBlock", "reached max limit maxOmegas", 1);
			}
			else if ( !(strcmp(currBlock, "HV")) ) {
				if ( newLag ) {
					j = model->p++;
					if ( j == 0 )
						if ( !(model->AR = (sparseMatrix **) arr_alloc(maxP, sparseMatrix *)) )
							errMsg("allocation", "readLinTrans", "AR", 0);
					if ( !(model->AR[j] = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))) )
						errMsg("allocation", "readLinTrans", "AR[n]", 0);
					if ( !(model->AR[j]->row = (iVector) arr_alloc(maxMatcnt, int)) )
						errMsg("allocation", "readLinTrans" ,"AR[n] rows", 0);
					if ( !(model->AR[j]->col = (iVector) arr_alloc(maxMatcnt, int)) )
						errMsg("allocation", "readLinTrans" ,"AR[n] columns", 0);
					if ( !(model->AR[j]->val = (dVector) arr_alloc(maxMatcnt, double)) )
						errMsg("allocation", "readLinTrans" ,"AR[n] coefficients", 0);
					model->AR[j]->cnt = 0;
				}
				row = 0;
				while ( (strcmp((*rvCols)[row], fields[0])) || (strcmp((*rvRows)[row], fields[1])) )
					row++;
				if ( row >= offset) {
					row -= model->M;
					period = 0;
					while ( period < numPeriods) {
						if( row < periodBeg[period] )
							break;
						period++;
					}
					row -= periodBeg[period-1];
				}
				model->AR[j]->row[model->AR[j]->cnt] = row;
				model->AR[j]->col[model->AR[j]->cnt] = col;
				model->AR[j]->val[model->AR[j]->cnt] = str2float(fields[2]);
				model->AR[j]->cnt++;
			}
			else if ( !(strcmp(currBlock, "LV")) || !(strcmp(currBlock, "RV")) ) {
				if ( newLag ) {
					j = model->q++;
					if ( j == 0 ) {
						if ( !(model->MA = (sparseMatrix **) arr_alloc(maxQ, sparseMatrix *)) )
							errMsg("allocation", "readLinTrans", "MA", 0);
					}
					if ( !(model->MA[j] = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))) )
						errMsg("allocation", "readLinTrans", "AR[n]", 0);
					if ( !(model->MA[j]->row = (iVector) arr_alloc(maxMatcnt, int)) )
						errMsg("allocation", "readLinTrans" ,"AR[n] rows", 0);
					if ( !(model->MA[j]->col = (iVector) arr_alloc(maxMatcnt, int)) )
						errMsg("allocation", "readLinTrans" ,"AR[n] columns", 0);
					if ( !(model->MA[j]->val = (dVector) arr_alloc(maxMatcnt, double)) )
						errMsg("allocation", "readLinTrans" ,"AR[n] coefficients", 0);
					model->MA[j]->cnt = 0;
				}

				row = 0;
				while ( (strcmp((*rvCols)[row], fields[0])) || (strcmp((*rvRows)[row], fields[1])) )
					row++;

				if ( row >= offset) {
					row -= model->M;
					period = 0;
					while ( period < numPeriods) {
						if( row < periodBeg[period] )
							break;
						period++;
					}
					row -= periodBeg[period-1];
				}

				model->MA[j]->row[model->MA[j]->cnt] = row;
				model->MA[j]->col[model->MA[j]->cnt] = col;
				model->MA[j]->val[model->MA[j]->cnt] = str2float(fields[2]);
				model->MA[j]->cnt++;
			}
			else {
				errMsg("read", "readLinTrans", "unknown block type encountered", 0);
				return 1;
			}
		}
	}

	/* Reallocate memory to exact values */
	for ( j = 0; j < model->p; j++ ) {
		model->AR[j]->col = (iVector) mem_realloc(model->AR[j]->col, model->AR[j]->cnt*sizeof(int));
		model->AR[j]->row = (iVector) mem_realloc(model->AR[j]->row, model->AR[j]->cnt*sizeof(int));
		model->AR[j]->val = (dVector) mem_realloc(model->AR[j]->val, model->AR[j]->cnt*sizeof(double));
	}
	model->AR = (sparseMatrix **) mem_realloc(model->AR, model->p*sizeof(sparseMatrix *));
	for ( j = 0; j < model->q; j++ ) {
		model->MA[j]->col = (iVector) mem_realloc(model->MA[j]->col, model->MA[j]->cnt*sizeof(int));
		model->MA[j]->row = (iVector) mem_realloc(model->MA[j]->row, model->MA[j]->cnt*sizeof(int));
		model->MA[j]->val = (dVector) mem_realloc(model->MA[j]->val, model->MA[j]->cnt*sizeof(double));
	}
	model->MA = (sparseMatrix **) mem_realloc(model->MA, model->q*sizeof(sparseMatrix *));

	/* Allocate trend as the mean value when available, else initial values as the mean values */
	if ( model->eta == NULL ) {
		for ( j = 0; j < stoc->numOmega; j++ )
			stoc->mean[j] = stoc->vals[0][j];
	}

	stoc->mod = model;
	mem_free(periodBeg);
	return 0;
}//END readLinTrans()

int readScenarios(FILE *fptr, cString *fields, oneProblem *orig, timeType *tim, int maxOmegas, int maxVals, stocType *stoc) {
	cString	*rvRows, *rvCols, *scenName;
	char 	strType;
	int  	n, r, c, numFields, maxScenarios = 100, numScen = 0, parentIdx;

	/* allocate memory to hold the names of random variable */
	if ( !(rvRows = (cString *) arr_alloc(maxOmegas, cString)) )
		errMsg("allocation", "readScenarios", "rvNames", 0);
	if ( !(rvCols = (cString *) arr_alloc(maxOmegas, cString)) )
		errMsg("allocation", "readScenarios", "rvNames", 0);
	if ( !(scenName = (cString *) arr_alloc(maxScenarios, cString)) )
		errMsg("allocation", "readScenarios", "scenNames", 0);

	stoc->vals     = (dVector *) arr_alloc(maxOmegas, dVector);
	stoc->numVals  = (iVector) arr_alloc(maxOmegas, int);
	stoc->probs    = (dVector *) arr_alloc(1, dVector);
	stoc->probs[0] = (dVector) arr_alloc(maxScenarios, double);

	if ( !(strcmp(fields[1], "DISCRETE")) ) {
		/* store the type of stochastic process encountered */
		sprintf(stoc->type, "SCENARIOS_DISCRETE");
		while (true) {
			NEXT_LINE: getLine(&fptr, fields, &strType, &numFields);
			if (strType != 'f')
				break; 										//Encountered ENDATA
			if ( !(strcmp(fields[0], "SC")) ) {
				/* New scenario encountered */
				if (!(strcmp(fields[2], "ROOT")) ) {
					/* The current scenario is the root scenario. Need to copy the names of random variable rows and columns. */
					if ( !(scenName[numScen] = (cString) arr_alloc(NAMESIZE, char)) )
						errMsg("allocation", "readScenarios", "scenNames[n]", 0);
					strcpy(scenName[numScen], fields[1]);
					parentIdx = -1;
				}
				else {
					/* a non-root scenario encountered */
					if ( !(scenName[numScen] = (cString) arr_alloc(NAMESIZE, char)) )
						errMsg("allocation", "readScenarios", "scenNames[n]", 0);
					strcpy(scenName[numScen], fields[1]);
					parentIdx = 0;
					while (parentIdx < numScen) {
						if ( !(strcmp(scenName[parentIdx], fields[2])) )
							break;
						parentIdx++;
					}
				}
				if ( parentIdx >= 0 )
					for ( n = 0; n < stoc->numOmega; n++ )
						stoc->vals[n][numScen] = stoc->vals[n][parentIdx];
				stoc->probs[0][numScen++] = str2float(fields[3]);
			}
			else {
				n = 0;
				while (n < stoc->numOmega ) {
					if ( !(strcmp(fields[0], rvCols[n])) && !(strcmp(fields[1], rvRows[n])) )
						break;
					n++;
				}
				if ( n == stoc->numOmega ) {
					/* New omega encountered. Note the row and column information. */
					/* First, begin by identifying row and column coordinates in the problem */
					if ( !(strcmp(fields[0], "RHS")) )
						c = -1;
					else {
						c = 0;
						while ( c < orig->mac ){
							if ( !(strcmp(fields[0], orig->cname[c])) )
								break;
							c++;
						}
					}
					if ( c == orig->mac ) {
						errMsg("read", "readIndep", "unknown column name in the stoch file", 0);
						return 1;
					}
					else if ( c < tim->col[1] && c != -1)
						/* variable belongs to the first stage, hence ignore */
						goto NEXT_LINE;

					if ( !(strcmp(fields[1], orig->objname)) )
						r = -1;
					else {
						r = 0;
						while (r < orig->mar ) {
							if ( !(strcmp(fields[1], orig->rname[r])) )
								break;
							r++;
						}
					}
					if ( r == orig->mar ) {
						errMsg("read", "readIndep", "unknown row name in the stoch file", 0);
						return 1;
					}
					else if ( r < tim->row[1] && r != -1)
						/* variable belongs to the first stage, hence ignore */
						goto NEXT_LINE;

					if ( !(rvRows[stoc->numOmega] = (cString) arr_alloc(NAMESIZE, char)) )
						errMsg("allocation", "readIndep", "rvNames[n]", 0);
					if ( !(rvCols[stoc->numOmega] = (cString) arr_alloc(NAMESIZE, char)) )
						errMsg("allocation", "readIndep", "rvNames[n]", 0);
					if ( !(stoc->vals[stoc->numOmega] = (dVector) arr_alloc(maxVals, double)) )
						errMsg("allocation", "readIndep","omega.vals[n]", 0);

					strcpy(rvCols[stoc->numOmega], fields[0]);
					strcpy(rvRows[stoc->numOmega], fields[1]);
					stoc->numVals[stoc->numOmega] = 0;

					stoc->col[stoc->numOmega] = c;
					stoc->row[stoc->numOmega] = r;
					stoc->numOmega++;
				}
				stoc->vals[n][numScen-1] = str2float(fields[2]);
				stoc->numVals[n]++;
			}
		}
	}
	else {
		errMsg("read", "readScenarios", "unknown keyword in stoch file", 0);
		return 1;
	}

	/* free up memory from temporary variables */
	for ( n = 0; n < stoc->numOmega; n++) {
		if ( rvRows[n] ) mem_free(rvRows[n]);
		if ( rvCols[n] ) mem_free(rvCols[n]);
	}
	mem_free(rvRows); mem_free(rvCols);
	for ( n = 0; n < numScen; n++)
		if (scenName[n]) mem_free(scenName[n]);
	mem_free(scenName);

	/* compute the mean value for each random variable */
	for (n = 0; n < stoc->numOmega; n++ ) {
		for ( c = 0; c < numScen; c++ )
			stoc->mean[n] += stoc->probs[0][c]*stoc->vals[n][c];
	}

	/* reallocate memory to elements of stocType based on the exact sizes */
	stoc->col 		= (iVector) mem_realloc(stoc->col, stoc->numOmega*sizeof(int));
	stoc->row 		= (iVector) mem_realloc(stoc->row, stoc->numOmega*sizeof(int));
	stoc->mean      = (dVector) mem_realloc(stoc->mean, stoc->numOmega*sizeof(double));
	stoc->numVals	= (iVector) mem_realloc(stoc->numVals, stoc->numOmega*sizeof(int));
	for ( n = 0; n < stoc->numOmega; n++ )
		stoc->vals[n] = (dVector) mem_realloc(stoc->vals[n], numScen*sizeof(double));
	stoc->vals 		= (dVector *) mem_realloc(stoc->vals, stoc->numOmega*sizeof(dVector));
	stoc->probs[0] 	= (dVector) mem_realloc(stoc->probs[0], numScen*sizeof(double));
	stoc->probs 	= (dVector *) mem_realloc(stoc->probs, 1*sizeof(dVector));
	mem_free(stoc->groupBeg); stoc->groupBeg = NULL;
	mem_free(stoc->numPerGroup);
	stoc->numGroups++;
	stoc->numPerGroup = NULL;

	return 0;
}//END readScenarios()

void freeOneProblem(oneProblem *p) {

	if(p){
		if ( p->lp ) freeProblem(p->lp);
		if(p->objname) mem_free(p->objname);
		if(p->name) mem_free(p->name);
		if(p->objx) mem_free(p->objx);
		if(p->rhsx) mem_free(p->rhsx);
		if(p->senx) mem_free(p->senx);
		if(p->bdl) mem_free(p->bdl);
		if(p->bdu) mem_free(p->bdu);
		if(p->ctype) mem_free(p->ctype);
		if(p->matbeg) mem_free(p->matbeg);
		if(p->matval) mem_free(p->matval);
		if(p->matind) mem_free(p->matind);
		if(p->matcnt) mem_free(p->matcnt);
		if(p->cname) mem_free(p->cname);
		if(p->rname) mem_free(p->rname);
		if(p->cstore) mem_free(p->cstore);
		if(p->rstore) mem_free(p->rstore);
		mem_free(p);
	}

}//END freeOneProblem()

void freeTimeType(timeType *tim) {
	int n;

	if(tim){
//		if (tim->probName) mem_free(tim->probName);
		if (tim->colStg) mem_free(tim->colStg);
		if (tim->rowStg) mem_free(tim->rowStg);
		if (tim->col) mem_free(tim->col);
		if (tim->row) mem_free(tim->row);
		if (tim->stgNames) {
			for (n = 0; n < tim->numStages; n++ )
				if (tim->stgNames[n]) mem_free(tim->stgNames[n]);
			mem_free(tim->stgNames);
		}
		mem_free(tim);
	}
}//END freeTime()

void freeStocType(stocType *stoc) {
	int n;

	if ( stoc ) {
		if ( stoc->col ) mem_free(stoc->col);
		if ( stoc->row ) mem_free(stoc->row);
		if ( stoc->mean ) mem_free(stoc->mean);
		if ( stoc->numVals ) mem_free(stoc->numVals);
		if ( stoc->vals) {
			if ( !(strcmp(stoc->type, "INDEP_NORMAL")) || !(strcmp(stoc->type, "LINTRAN")))
				mem_free(stoc->vals[0]);
			else {
				for ( n = 0; n < stoc->numOmega; n++ )
					if ( stoc->vals[n] ) mem_free(stoc->vals[n]);
			}
			mem_free(stoc->vals);
		}
		if ( stoc->probs ) {
			if ( !(strcmp(stoc->type, "SCENARIOS_DISCRETE")) ) {
				if (stoc->probs[0]) mem_free(stoc->probs[0]);
			}
			else {
				for ( n = 0; n < maximum(stoc->numOmega, stoc->numGroups); n++ )
					if ( stoc->probs[n] ) mem_free(stoc->probs[n]);
			}
			mem_free(stoc->probs);
		}
		if ( stoc->groupBeg) mem_free(stoc->groupBeg);
		if ( stoc->numPerGroup) mem_free(stoc->numPerGroup);
		if ( stoc->type) mem_free(stoc->type);
		if ( stoc->mod) freeStatModel(stoc->mod);
		mem_free(stoc);
	}

}//END freeStocType()

void freeStatModel(statModel *model) {
	int n;

	if ( model->AR ) {
		for ( n = 0; n < model->p; n++ )
			if (model->AR[n]) freeSparseMatrix(model->AR[n]);
		mem_free(model->AR);
	}
	if ( model->MA ) {
		for ( n = 0; n < model->q; n++ )
			if (model->MA[n]) freeSparseMatrix(model->MA[n]);
		mem_free(model->MA);
	}
	if ( model->eta ) {
		for ( n = 0; n < model->N; n++ )
			if ( model->eta[n]) mem_free(model->eta[n]);
		mem_free(model->eta);
	}
	if ( model->sigma ) {
		for ( n = 0; n < model->N; n++ )
			if ( model->sigma[n]) mem_free(model->sigma[n]);
		mem_free(model->sigma);
	}
	if ( model->muEps ) mem_free(model->muEps);
	if ( model->cvEps ) freeSparseMatrix(model->cvEps);

	mem_free(model);

}//END freeStatModel()
