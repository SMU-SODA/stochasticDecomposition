/*
 * setup.c
 *
 *  Created on: Mar 27, 2017
 *      Author: gjharsha
 */

#include "sddp.h"

extern configType config;

int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType ***cell) {
	vector	meanSol, lb;
	int 	t;

	/* setup mean value problem which will act as reference for all future computations */
	meanSol = meanProblem(orig, stoc);
	if ( meanSol == NULL ) {
		errMsg("setup", "setupAlgo", "failed to setup and solve mean value problem", 0);
		return 1;
	}

	/* calculate lower bounds for each stage */
	lb = calcLowerBound(orig, tim, stoc);
	if ( lb == NULL ) {
		errMsg("setup", "setupAlgo", "failed to compute lower bounds on stage problem", 0);
		return 1;
	}

	/* decompose the problem into master and subproblem */
	(*prob) = newProb(orig, stoc, tim, lb, config.TOLERANCE);
	if ( (*prob) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to update probType with elements specific to algorithm", 0);
		return 1;
	}

#ifdef DECOMPOSE_CHECK
	printDecomposeSummary(tim, (*prob));
#endif

	/* ensure that we have a linear programs at all stages */
	t = 0;
	while ( t < tim->numStages ) {
		if ( (*prob)[t++]->sp->type  != PROB_LP )
			printf("Warning :: Stage-%d problem is a mixed-integer program. Solving its linear relaxation.\n", t);
	}

	/* create the cells which will be used in the algorithms */
	(*cell) = newCell(stoc, (*prob), lb, tim->numStages);
	if ( (*cell) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to create the necessary cell structure", 0);
		return 1;
	}

	mem_free(meanSol); mem_free(lb);

	return 0;
}//END setupAlgo()

cellType **newCell(stocType *stoc, probType **prob, vector lb, int T) {
	cellType 	**cell;
	int			t, i, j, cnt, rOffset, cOffset;
	char		*q;

	if ( !(cell = (cellType **) arr_alloc(T, cellType *)) )
		errMsg("allocation", "newCell", "cell", 0);

	for ( t = 0; t < T; t++ ) {
		if ( !(cell[t] = (cellType *) mem_malloc(sizeof(cellType))) )
			errMsg("allocation", "newCell", "cell[t]", 0);
		cell[t]->k = 0;
		cell[t]->lb = lb[t];
		if ( DBL_ABS(lb[t]) > config.TOLERANCE )
			cell[t]->lbType = NONTRIVIAL;
		else
			cell[t]->lbType = TRIVIAL;

		/* oneProblem type */
		if ( !(cell[t]->sp = (oneProblem *) mem_malloc(sizeof(oneProblem))) )
			errMsg("allocation", "newCell", "cell[t]->sp", 0);

		/* assign values from probType */
		cell[t]->sp->matsz = prob[t]->sp->matsz;
		cell[t]->sp->macsz = prob[t]->sp->mac;
		cell[t]->sp->mac = prob[t]->sp->mac;
		cell[t]->sp->cstorsz = prob[t]->sp->cstorsz;
		cell[t]->sp->marsz = prob[t]->sp->mar;
		cell[t]->sp->mar = prob[t]->sp->mar;
		cell[t]->sp->rstorsz = prob[t]->sp->rstorsz;
		cell[t]->sp->objsen = prob[t]->sp->objsen;
		cell[t]->sp->type = PROB_LP;
		if ( t != T-1 ) {
			cell[t]->sp->macsz++;
			cell[t]->sp->mac++;
			cell[t]->sp->cstorsz += NAMESIZE;
		}

		if ( !(cell[t]->sp->name = (string) arr_alloc(NAMESIZE, char)) )
			errMsg("allocation", "newCell", "cell[t]->sp->name", 0);
		if ( !(cell[t]->sp->objname = (string) arr_alloc(NAMESIZE, char)) )
			errMsg("allocation", "newCell", "cell[t]->sp->objname", 0);
		if ( !(cell[t]->sp->objx = (vector) arr_alloc(cell[t]->sp->mac, double)) )
			errMsg("allocation", "newCell", "cell[t]->sp->objx", 0);
		if ( !(cell[t]->sp->bdl = (vector) arr_alloc(cell[t]->sp->mac, double)) )
			errMsg("allocation", "newCell", "cell[t]->sp->bdl", 0);
		if ( !(cell[t]->sp->bdu = (vector) arr_alloc(cell[t]->sp->mac, double)) )
			errMsg("allocation", "newCell", "cell[t]->sp->bdu", 0);
		if ( !(cell[t]->sp->matbeg = (intvec) arr_alloc(cell[t]->sp->mac, int)) )
			errMsg("allocation", "newCell", "cell[t]->sp->matbeg", 0);
		if ( !(cell[t]->sp->matcnt = (intvec) arr_alloc(cell[t]->sp->mac, int)) )
			errMsg("allocation", "newCell", "cell[t]->sp->matcnt", 0);
		if ( !(cell[t]->sp->rhsx = (vector) arr_alloc(cell[t]->sp->mar, double)) )
			errMsg("allocation", "newCell", "cell[t]->sp->rhsx", 0);
		if ( !(cell[t]->sp->senx = (string) arr_alloc(cell[t]->sp->mar, char)) )
			errMsg("allocation", "newCell", "cell[t]->sp->senx", 0);
		if (!(cell[t]->sp->cname = (string *) arr_alloc(cell[t]->sp->macsz, string)))
			errMsg("Allocation", "newCell", "cell[t]->sp->cname",0);
		if (!(cell[t]->sp->cstore = (string) arr_alloc(cell[t]->sp->cstorsz, char)))
			errMsg("Allocation", "newCell", "cell[t]->sp->cstore", 0);
		if (!(cell[t]->sp->rname = (string *)arr_alloc(cell[t]->sp->marsz, string)))
			errMsg("Allocation", "newCell", "cell[t]->sp->rname", 0);
		if (!(cell[t]->sp->rstore = (string) arr_alloc(cell[t]->sp->rstorsz, char)))
			errMsg("Allocation", "newCell", "cell[t]->sp->rstore", 0);
		if (!(cell[t]->sp->matval = (vector) arr_alloc(cell[t]->sp->matsz, double)))
			errMsg("Allocation", "newCell", "cell[t]->sp->matval", 0);
		if (!(cell[t]->sp->matind = (intvec) arr_alloc(cell[t]->sp->matsz, int)))
			errMsg("Allocation", "newCell", "cell[t]->sp->matind", 0);
		cell[t]->sp->ctype = NULL;

		strcpy(cell[t]->sp->name, prob[t]->sp->name);
		strcpy(cell[t]->sp->objname, prob[t]->sp->objname);

		/* copy the cell problem's column and row names */
		i = 0;
		for (q = prob[t]->sp->cname[0]; q < prob[t]->sp->cname[0] + prob[t]->sp->cstorsz; q++)
			cell[t]->sp->cstore[i++] = *q;

		i = 0;
		for (q = prob[t]->sp->rname[0]; q < prob[t]->sp->rname[0] + prob[t]->sp->rstorsz; q++)
			cell[t]->sp->rstore[i++] = *q;

		/* Calculate difference in pointers for cell cell[t]-> row and column names */
		cOffset = cell[t]->sp->cstore - prob[t]->sp->cname[0];
		rOffset = cell[t]->sp->rstore - prob[t]->sp->rname[0];

		/* Copy the all column information from problem in probType */
		cnt = 0;
		for (i = 0; i < prob[t]->sp->mac; i++)	{
			cell[t]->sp->objx[i] = prob[t]->sp->objx[i];
			cell[t]->sp->bdu[i] = prob[t]->sp->bdu[i];
			cell[t]->sp->bdl[i] = prob[t]->sp->bdl[i];
			cell[t]->sp->cname[i] = prob[t]->sp->cname[i] + cOffset;
			cell[t]->sp->matbeg[i] = cnt;
			cell[t]->sp->matcnt[i] = prob[t]->sp->matcnt[i];
			for (j = prob[t]->sp->matbeg[i]; j < prob[t]->sp->matbeg[i] + prob[t]->sp->matcnt[i]; j++) {
				cell[t]->sp->matval[cnt] = prob[t]->sp->matval[j];
				cell[t]->sp->matind[cnt] = prob[t]->sp->matind[j];
				cnt++;
			}
		}

		/* Copy all information concerning rows of cell */
		for (i = 0; i < prob[t]->sp->mar; i++) {
			cell[t]->sp->rhsx[i] = prob[t]->sp->rhsx[i];
			cell[t]->sp->senx[i] = prob[t]->sp->senx[i];
			cell[t]->sp->rname[i] = prob[t]->sp->rname[i] + rOffset;
		}

		if ( t != T-1 ) {
			strcpy(cell[t]->sp->cstore + prob[t]->sp->cstorsz, "eta");
			cell[t]->sp->cname[prob[t]->sp->mac] = cell[t]->sp->cstore + prob[t]->sp->cstorsz;
			cell[t]->sp->objx[prob[t]->sp->mac] = 1.0;
			cell[t]->sp->bdu[prob[t]->sp->mac] = INFBOUND;
			cell[t]->sp->bdl[prob[t]->sp->mac] = cell[t]->lb;
			cell[t]->sp->matbeg[prob[t]->sp->mac] = cnt;
			cell[t]->sp->matcnt[prob[t]->sp->mac] = 0;
		}

		/* Solution part of the cell structure: primal and dual. */
		if ( t != T-1 ) {
			/* candidate solution */
			if ( !(cell[t]->candidU = (vector) arr_alloc(prob[t]->num->cols+1, double)) )
				errMsg("allocation", "newCell", "cell[t]->candidU", 0);
			/* TODO: incumbent solutions */
			cell[t]->candidEst = 0.0;

			/* cuts structure */
			if ( !(cell[t]->cuts = (cutsType *) mem_malloc(sizeof(cutsType))) )
				errMsg("allocation", "newCell", "cell[t]->cuts", 0);
			if ( !(cell[t]->cuts->vals = (oneCut **) arr_alloc(config.MAX_ITER, oneCut *)) )
				errMsg("allocation", "newCell", "cell[t]->cuts->vals", 0);
			cell[t]->cuts->cnt = 0;
		}
		else {
			if ( !(cell[t]->candidU = (vector) arr_alloc(prob[t]->num->cols+1, double)) )
				errMsg("allocation", "newCell", "cell[t]->candidU", 0);
			cell[t]->incumbU = NULL;
			cell[t]->cuts 	 = NULL;
		}

		if ( t > 0 ) {
			/* for non-root stages */
			if (!(cell[t]->pi = (vector) arr_alloc(prob[t]->num->rows+config.MAX_ITER+1, double)) )
				errMsg("allocation", "newCell", "cell[t]->pi", 0);
			if (!(cell[t]->rhs = (vector) arr_alloc(prob[t]->num->rows+1, double)) )
				errMsg("allocation", "newCell", "cell[t]->pi", 0);

			/* stochastic elements of the cell */
			cell[t]->omega = newOmega(t-1, stoc);
			if ( cell[t]->omega == NULL ) {
				errMsg("setup", "newCell", "failed to setup new omega structure", 0);
				return NULL;
			}

#if 0
			cell[t]->lambda = newLambda(cell[t]->omega->cnt);
			cell[t]->sigma 	= newSigma(cell[t]->omega->cnt, 0);
			cell[t]->delta = newDelta(cell[t]->omega->cnt);
#else if
			cell[t]->lambda = newLambda(config.MAX_ITER);
			cell[t]->sigma 	= newSigma(config.MAX_ITER, 0);
			cell[t]->delta = newDelta(config.MAX_ITER);
#endif

		}
		else {
			cell[t]->pi 	= NULL;
			cell[t]->rhs 	= NULL;
			cell[t]->omega	= NULL;
			cell[t]->lambda = NULL;
			cell[t]->sigma 	= NULL;
			cell[t]->delta  = NULL;
		}

		/* Load the cell problem onto solver: this problem will be used in forward pass */
		cell[t]->sp->lp = setupProblem(cell[t]->sp->name, cell[t]->sp->type, cell[t]->sp->mac, cell[t]->sp->mar,
				cell[t]->sp->objsen, cell[t]->sp->objx, cell[t]->sp->rhsx, cell[t]->sp->senx, cell[t]->sp->matbeg,
				cell[t]->sp->matcnt, cell[t]->sp->matind, cell[t]->sp->matval, cell[t]->sp->bdl, cell[t]->sp->bdu,
				NULL, cell[t]->sp->cname, cell[t]->sp->rname, cell[t]->sp->ctype);
		if ( cell[t]->sp->lp == NULL ) {
			errMsg("solver", "newCell", "failed to setup cell problem on solver",0);
			return NULL;
		}

#if CELL_SETUP
		char fname[NAMESIZE];
		sprintf(fname, "cell%d.lp", t);
		writeProblem(cell[t]->sp->lp, fname);
#endif
	}

	return cell;
}//END newCell()

void cleanupAlgo(probType **prob, cellType **cell, int T) {

	freeCellType(prob, cell, T);
	freeProbType(prob, T);

}//END cleanupAlgo()

void freeCellType (probType **prob, cellType **cell, int T) {
	int t;

	if (cell) {
		for ( t = 0; t < T; t++ ) {
			if ( cell[t] ) {
				if (cell[t]->sp) freeOneProblem(cell[t]->sp);
				if (cell[t]->candidU) mem_free(cell[t]->candidU);
				if (cell[t]->pi) mem_free(cell[t]->pi);
				if (cell[t]->rhs) mem_free(cell[t]->rhs);
				if (cell[t]->cuts) freeCutsType(cell[t]->cuts);
				if (cell[t]->lambda) freeLambdaType(cell[t]->lambda, TRUE);
				if (cell[t]->sigma) freeSigmaType(cell[t]->sigma, TRUE);
				if (cell[t]->delta) freeDeltaType(cell[t]->delta, cell[t]->omega->cnt, TRUE);
				if (cell[t]->omega) freeOmegaType(cell[t]->omega);
				mem_free (cell[t]);
			}
		}
		mem_free(cell);
	}

}//END freeCellType
