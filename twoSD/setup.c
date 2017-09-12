/*
 * setup.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "cell.h"

extern configType config;

int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell) {
	vector	meanSol = NULL, lb = NULL;
	int 	t;

	/* setup mean value problem which will act as reference for all future computations */
	meanSol = meanProblem(orig, stoc);
	if ( meanSol == NULL ) {
		errMsg("setup", "setupAlgo", "failed to setup and solve mean value problem", 0);
		goto TERMINATE;
	}

	/* calculate lower bounds for each stage */
	lb = calcLowerBound(orig, tim, stoc);
	if ( lb == NULL )  {
		errMsg("setup", "setupAlgo", "failed to compute lower bounds on stage problem", 0);
		goto TERMINATE;
	}

	/* decompose the problem into master and subproblem */
	(*prob) = newProb(orig, stoc, tim, lb, config.TOLERANCE);
	if ( (*prob) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to update probType with elements specific to algorithm", 0);
		goto TERMINATE;
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
	(*cell) = newCell(stoc, (*prob), meanSol);
	if ( (*cell) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to create the necessary cell structure", 0);
		goto TERMINATE;
	}

	mem_free(meanSol); mem_free(lb);
	return 0;
	TERMINATE:
	mem_free(meanSol); mem_free(lb);
	return 1;
}//END setupAlgo()

/* This function is used to create cells used in the algorithm */
cellType *newCell(stocType *stoc, probType **prob, vector xk) {
	cellType    *cell;
	int			length;

	/* Allocate memory to all cells used in the algorithm. */
	if (!(cell = (cellType *) mem_malloc(sizeof(cellType))) )
		errMsg("Memory allocation", "new_cell", "failed to allocate memory to cell",0);
	cell->master = cell->subprob = NULL;
	cell->candidX = cell->incumbX = NULL;
	cell->piM = cell->djM = NULL;
	cell->cuts = cell->fcuts = NULL;
	cell->lambda = NULL; cell->sigma = NULL; cell->delta = NULL; cell->omega = NULL; cell->basis = NULL;
	cell->pi_ratio = NULL;

	/* setup the master problem */
	cell->master = newMaster(prob[0]->sp, prob[0]->lb);
	if ( cell->master == NULL ) {
		errMsg("setup", "newCell", "failed to setup the master problem", 0);
		return NULL;
	}
	/* setup the subproblem */
	cell->subprob = newSubproblem(prob[1]->sp);

	/* -+-+-+-+-+-+-+-+-+-+-+ Allocating memory to other variables that belongs to cell +-+-+-+-+-+-+-+-+-+- */
	cell->k 	= 0;
	cell->LPcnt = 0;
	if (prob[0]->lb == 0)
		cell->lbType = TRIVIAL;
	else
		cell->lbType = NONTRIVIAL;
	cell->lb = prob[0]->lb;

	/* candidate solution and estimates */
	cell->candidX 	= duplicVector(xk, prob[0]->num->cols);
	cell->candidEst	= prob[0]->lb + vXvSparse(cell->candidX, prob[0]->dBar);
	cell->ub 		= INF;

	/* incumbent solution and estimates */
	if (config.MASTER_TYPE == PROB_QP) {
		cell->incumbX   = duplicVector(xk, prob[0]->num->cols);
		cell->incumbEst = cell->candidEst;
		cell->quadScalar= config.MIN_QUAD_SCALAR;     						/* The quadratic scalar, 'sigma'*/
		cell->iCutIdx   = 0;
		cell->iCutUpdt  = 0;
		cell->incumbChg = TRUE;
	}
	else {
		cell->incumbX   = NULL;
		cell->incumbEst = 0.0;
		cell->quadScalar= 0.0;
		cell->iCutIdx   = -1;
		cell->iCutUpdt  = -1;
		cell->incumbChg = FALSE;
	}
	cell->gamma 			= 0.0;
	cell->normDk_1 			= 0.0;
	cell->normDk 			= 0.0;

	/* lower bounding approximations held in cuts structure */
	cell->maxCuts = config.CUT_MULT * prob[0]->num->cols + 3;
	cell->cuts 	  = newCuts(cell->maxCuts);
	cell->fcuts   = NULL;

	/* solution parts of the cell */
	if ( !(cell->djM = (vector) arr_alloc(prob[0]->num->cols + 2, double)) )
		errMsg("allocation", "newMaster", "cell->di", 0);
	if ( !(cell->piM = (vector) arr_alloc(prob[0]->num->rows + cell->maxCuts + 1, double)) )
		errMsg("allocation", "newMaster", "cell->piM", 0);

	/* stochastic elements */
	if ( prob[1]->num->rvdOmCnt > 0 )
		length = prob[1]->num->rvdOmCnt*config.MAX_ITER + config.MAX_ITER / config.TAU + 1;
	else
		length = config.MAX_ITER + config.MAX_ITER / config.TAU + 1;
	cell->basis  = newBasisType(prob[1]->sp->senx, config.MAX_ITER, prob[1]->num->cols, prob[1]->num->rows, WORDLENGTH);
	cell->lambda = newLambda(length, 0, prob[1]->num->rvRowCnt);
	cell->sigma  = newSigma(length, prob[1]->num->cntCcols, 0);
	cell->delta  = newDelta(length);
	cell->omega  = newOmega(stoc->numOmega, config.MAX_ITER);

	cell->optFlag 			= FALSE;
	cell->dualStableFlag 	= FALSE;
	if ( !(cell->pi_ratio = (vector) arr_alloc(config.SCAN_LEN, double)) )
		errMsg("allocation", "newCell", "cell->pi_ratio", 0);

	cell->feasCnt 			= 0;
	cell->infeasIncumb 		= FALSE;

	/* construct the QP using the current incumbent */
	if ( config.MASTER_TYPE == PROB_QP ) {
		if ( constructQP(prob[0], cell, cell->incumbX, cell->quadScalar) ) {
			errMsg("setup", "newCell", "failed to change the right-hand side after incumbent change", 0);
			return NULL;
		}

		cell->incumbChg = FALSE;
#if defined(SETUP_CHECK)
		if ( writeProblem(cell->aster->lp, "newQPMaster.lp") ) {
			errMsg("write problem", "new_master", "failed to write master problem to file",0);
			return NULL;
		}
#endif
	}

	return cell;
}//END newCell()

int readConfig() {
	FILE 	*fptr;
	char	line[2*BLOCKSIZE], comment[2*BLOCKSIZE];
	int 	status;

	fptr = fopen("config.sd", "r");
	if ( fptr == NULL ) {
		errMsg("read", "readConfig", "failed to open configuration file", 0);
		return 1;
	}

	while ((status = (fscanf(fptr, "%s", line) != EOF))) {
		if (!(strcmp(line, "RUN_SEED")))
			fscanf(fptr, "%lld", &config.RUN_SEED);
		else if (!(strcmp(line, "TOLERANCE")))
			fscanf(fptr, "%lf", &config.TOLERANCE);
		else if (!(strcmp(line, "MIN_ITER")))
			fscanf(fptr, "%d", &config.MIN_ITER);
		else if (!(strcmp(line, "MAX_ITER")))
			fscanf(fptr, "%d", &config.MAX_ITER);
		else if (!(strcmp(line, "MASTER_TYPE")))
			fscanf(fptr, "%d", &config.MASTER_TYPE);
		else if (!(strcmp(line, "TAU")))
			fscanf(fptr, "%d", &config.TAU);
		else if (!(strcmp(line, "MIN_QUAD_SCALAR")))
			fscanf(fptr, "%lf", &config.MIN_QUAD_SCALAR);
		else if (!(strcmp(line, "EPSILON")))
			fscanf(fptr, "%lf", &config.EPSILON);
		else if (!(strcmp(line, "EVAL_FLAG")))
			fscanf(fptr, "%d", &config.EVAL_FLAG);
		else if (!(strcmp(line, "EVAL_SEED")))
			fscanf(fptr, "%lld", &config.EVAL_SEED);
		else if (!(strcmp(line, "EVAL_MIN_ITER")))
			fscanf(fptr, "%d", &config.EVAL_MIN_ITER);
		else if (!(strcmp(line, "EVAL_ERROR")))
			fscanf(fptr, "%lf", &config.EVAL_ERROR);

		else if (!(strcmp(line, "CUT_MULT")))
			fscanf(fptr, "%d", &config.CUT_MULT);
		else if (!(strcmp(line, "MAX_QUAD_SCALAR")))
			fscanf(fptr, "%lf", &config.MAX_QUAD_SCALAR);
		else if (!(strcmp(line, "R1")))
			fscanf(fptr, "%lf", &config.R1);
		else if (!(strcmp(line, "R2")))
			fscanf(fptr, "%lf", &config.R2);
		else if (!(strcmp(line, "R3")))
			fscanf(fptr, "%lf", &config.R3);
		else if (!(strcmp(line, "PI_EVAL_START")))
			fscanf(fptr, "%d", &config.PI_EVAL_START);
		else if (!(strcmp(line, "PI_CYCLE")))
			fscanf(fptr, "%d", &config.PI_CYCLE);
		else if (!(strcmp(line, "BOOTSTRAP_REP")))
			fscanf(fptr, "%d", &config.BOOTSTRAP_REP);
		else if (!(strcmp(line, "PERCENT_PASS")))
			fscanf(fptr, "%lf", &config.PERCENT_PASS);
		else if (!(strcmp(line, "SCAN_LEN")))
			fscanf(fptr, "%d", &config.SCAN_LEN);
		else if (!(strcmp(line, "PRE_EPSILON")))
			fscanf(fptr, "%lf", &config.PRE_EPSILON);

		else if (!(strcmp(line, "MAX_OBS")))
			fscanf(fptr, "%d", &config.MAX_OBS);
		else if (!(strcmp(line, "SUBPROB_SAMPLE_PCT")))
			fscanf(fptr, "%lf", &config.SUBPROB_SAMPLE_PCT);
		else if (!(strcmp(line, "SUBPROB_SAMPLE_SEED")))
			fscanf(fptr, "%lld", &config.SUBPROB_SAMPLE_SEED);
		else if (!(strcmp(line, "SAA")))
			fscanf(fptr, "%d", &config.SAA);
		else if (!strcmp(line, "//"))
			fgets(comment, 2*BLOCKSIZE, fptr);
		else {
			printf ("%s\n", line);
			errMsg("read", "readConfig", "unrecognized parameter in configuration file", 1);
		}
	}

	fclose(fptr);

	return 0;
}//END readConfig()

void freeCellType(cellType *cell) {

	if ( cell ) {
		if (cell->master) freeOneProblem(cell->master);
		if (cell->candidX) mem_free(cell->candidX);
		if (cell->incumbX) mem_free(cell->incumbX);
		if (cell->piM) mem_free(cell->piM);
		if (cell->djM) mem_free(cell->djM);
		if (cell->cuts) freeCutsType(cell->cuts);
		if (cell->fcuts) freeCutsType(cell->fcuts);
		if (cell->basis) freeBasisType(cell->basis);
		if (cell->delta) freeDeltaType(cell->delta, cell->lambda->cnt, cell->omega->cnt);
		if (cell->omega) freeOmegaType(cell->omega);
		if (cell->lambda) freeLambdaType(cell->lambda);
		if (cell->sigma) freeSigmaType(cell->sigma);
		if (cell->pi_ratio) mem_free(cell->pi_ratio);
		mem_free(cell);
	}

}//END freeCellType()

/* This subroutine initializes the master problem by copying information from the decomposed prob[0](type: oneProblem) and adding a column for
 * theta for modified benders decomposition. */
oneProblem *newMaster(oneProblem *orig, double lb) {
	oneProblem 	*master;
	int         r, i, j, idx, cnt;
	long        colOffset, rowOffset;
	char        *q;

	if (!(master = (oneProblem *) mem_malloc (sizeof(oneProblem))))
		errMsg("Memory allocation", "new_master", "Faile to allocate memory to mcell->sp", 0);

	/* -+-+-+-+-+-+-+-+-+-+-+-+-+-+- Allocating memory to master -+-+-+-+-+-+-+-+-+-+-+-+-+-+- */
	master->type 	= config.MASTER_TYPE;                  	/* type of problem: LP, QP, MIP or MIQP */
	master->objsen 	= orig->objsen;                 	/* sense of the objective: 1 for minimization and -1 for maximization */
	master->mar 	= orig->mar;                       	/* number of rows */
	master->numInt 	= orig->numInt;                 	/* number of integer variables in the problem  */
	master->numnz 	= orig->numnz;                   	/* number of non-zero elements in constraint matrix */
	master->matsz 	= orig->matsz;                   	/* extended matrix size */
	master->marsz 	= orig->marsz;                   	/* extended row size */
	master->rstorsz = orig->rstorsz;               		/* memory size for storing row names */
	master->mac 	= orig->mac+1;           			/* number of columns + etas */
	master->macsz 	= orig->macsz + 1;       			/* extended column size */
	master->cstorsz 	= orig->cstorsz + NAMESIZE;    	/* memory size for storing column names */

	/* Allocate memory to the information whose type is string */
	if (!(master->name = (string) arr_alloc(NAMESIZE, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->name",0);
	if (!(master->senx = (string) arr_alloc(master->marsz,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->senx",0);
	if (!(master->ctype = (string) arr_alloc(master->macsz,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->ctype",0);
	if (!(master->objname = (string) arr_alloc(NAMESIZE,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->objname",0);
	if (!(master->rname = (string *) arr_alloc(master->marsz,string)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->rname",0);
	if (!(master->rstore = (string) arr_alloc(master->rstorsz, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->rstore",0);
	if (!(master->cname = (string*) arr_alloc(master->macsz,string)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->cname",0);
	if (!(master->cstore = (string) arr_alloc(master->cstorsz, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->cstore",0);

	/* Allocate memory to the information whose type is vector */
	if (!(master->objx = (vector) arr_alloc(master->macsz, double)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->objx",0);
	if (!(master->rhsx = (vector) arr_alloc(master->marsz, double)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->rhsx",0);
	if (!(master->matval = (vector) arr_alloc(master->matsz, double)))
		errMsg("allocation", "new_master", "master->matval",0);
	if (!(master->bdl = (vector) arr_alloc(master->macsz, double)))
		errMsg("allocation", "new_master", "master->bdl",0);
	if (!(master->bdu = (vector) arr_alloc(master->macsz, double)))
		errMsg("allocation", "new_master", "master->bdu",0);

	/* Allocate memory to the information whose type is intvec */
	if (!(master->matbeg = (intvec) arr_alloc(master->macsz, int)))
		errMsg("allocation", "new_master", "master->matbeg",0);
	if (!(master->matcnt = (intvec) arr_alloc(master->macsz, int)))
		errMsg("allocation", "new_master", "master->matcnt",0);
	if (!(master->matind = (intvec) arr_alloc(master->matsz, int)))
		errMsg("allocation", "new_master", "master->matind",0);

	strcpy(master->name, orig->name);           /* Copy problem name */
	strcpy(master->objname, orig->objname);     /* Copy objective name */

	/* Copy problem's column and row names */
	i = 0;
	for (q = orig->cname[0]; q < orig->cname[0] + orig->cstorsz; q++)
		master->cstore[i++] = *q;

	i = 0;
	for (q = orig->rname[0]; q < orig->rname[0] + orig->rstorsz; q++)
		master->rstore[i++] = *q;

	/* Calculate difference in pointers for master/copy row and column names */
	colOffset = master->cstore - orig->cname[0];
	rowOffset = master->rstore - orig->rname[0];

	/* Copy the all column information from the original master problem */
	cnt = 0;
	for (j = 0; j < orig->mac; j++) {
		/* Copy objective function coefficients */
		master->objx[j] = orig->objx[j];
		/* Copy the decision variable type */
		master->ctype[j] = orig->ctype[j];
		/* Copy the upper bound and lower bound */
		master->bdu[j] = orig->bdu[j];
		master->bdl[j] = orig->bdl[j];
		/* Copy column names, offset by length */
		master->cname[j] = orig->cname[j] + colOffset;
		/* Copy the master sparse matrix beginning position of each column */
		master->matbeg[j] = cnt;
		/* Copy the sparse matrix non-zero element count */
		master->matcnt[j] = orig->matcnt[j];
		master->ctype[j] = orig->ctype[j];
		/* Loop through all non-zero elements in this column */
		for (idx = orig->matbeg[j]; idx < orig->matbeg[j] + orig->matcnt[j]; idx++) {
			/* Copy the non-zero coefficient */
			master->matval[cnt] = orig->matval[idx];
			/* Copy the row entry of the non-zero elements */
			master->matind[cnt] = orig->matind[idx];
			cnt++;
		}
	}

	/* Copy all information concerning rows of master */
	for (r = 0; r < orig->mar; r++) {
		/* Copy the right hand side value */
		master->rhsx[r] = orig->rhsx[r];
		/* Copy the constraint sense */
		master->senx[r] = orig->senx[r];
		/* Copy row names, offset by length */
		master->rname[r] = orig->rname[r] + rowOffset;
	}

	/* Initialize information for the extra column in the new master. */
	colOffset = orig->cstorsz;
	strcpy(master->cstore + orig->cstorsz, "eta");
	master->cname[orig->mac] = master->cstore + colOffset;
	master->objx[orig->mac] = 1.0;			// orig->mac is the last column in the original master
	master->ctype[orig->mac] = 'C';
	master->bdu[orig->mac] = INFBOUND;
	master->bdl[orig->mac] = lb;
	master->matbeg[orig->mac] = orig->numnz;	// Beginning point in matval/matind in eta columns. every eta column begins at the same address
	master->matcnt[orig->mac] = 0;               // Only optimality cuts has eta

	/* Load the copy into CPLEX */
	master->lp = setupProblem(master->name, master->type, master->mac, master->mar, master->objsen, master->objx, master->rhsx, master->senx, master->matbeg, master->matcnt,master->matind, master->matval, master->bdl, master->bdu, NULL, master->cname, master->rname, master->ctype);
	if ( master->lp == NULL ) {
		errMsg("Problem Setup", "new_master", "failed to setup master problem in the solver",0);
		return NULL;
	}

#if defined(SETUP_CHECK)
	if ( writeProblem(master->lp, "newMaster.lp") ) {
		errMsg("solver", "newMaster", "failed to write master problem to file", 0);
		return NULL;
	}
#endif

	return master;

}//END newMaster()

int constructQP(probType *prob, cellType *cell, vector incumbX, double quadScalar) {
	int status;

	status = changeQPproximal(cell->master->lp, prob->num->cols, quadScalar);
	if ( status ) {
		errMsg("algorithm", "algoIntSD", "failed to change the proximal term", 0);
		return 1;
	}
	status = changeQPrhs(prob, cell, incumbX);
	if ( status ) {
		errMsg("algorithm", "algoIntSD", "failed to change the right-hand side to convert the problem into QP", 0);
		return 1;
	}
	status = changeQPbds(cell->master->lp, prob->num->cols, prob->sp->bdl, prob->sp->bdu, incumbX);
	if ( status ) {
		errMsg("algorithm", "algoIntSD", "failed to change the bounds to convert the problem into QP", 0);
		return 1;
	}

	return status;
}//END constructQP()

/* Construct the Q diagonal matrix and copy it for quadratic problem. */
int changeQPproximal(LPptr lp, int numCols, double sigma) {
	int    n;
	vector qsepvec;

	if (!(qsepvec = arr_alloc(numCols+1, double)))
		errMsg("Allocation", "constructQP", "qsepvec",0);

	/* Construct Q matrix, which is simply a diagonal matrix. */
	for (n = 0; n < numCols; n++)
		qsepvec[n] = 0.5 * sigma;
	qsepvec[n] = 0.0;

	/* Now copy the Q matrix for QP problem. */
	if ( copyQPseparable(lp, qsepvec) ) {
		errMsg("solver", "constructQP", "failed to copy Q matrix", 0);
		return 1;
	}

	mem_free(qsepvec);
	return 0;
}//END changeQPproximal()

/* In the regularized QP method, we need to change the rhs of x to d. The
 * 		 A * x 			= b
 * 		 eta + beta * x >= alpha
 * Since x = xbar + d, the corresponding changes will therefore be:
 * 		 A * d = b - A * xbar
 * 		 eta + beta * d >= alpha - beta * xbar
 * But as long as the incumbent sulotion does not change, b - A * xbar and alpha - beta * xbar (for the existing cuts) won't change. So we only need
 * to change it when the incumbent changes.
 *
 * On the other hand, in each iteration, a new cut will be added (and/or some cuts may be dropped) and therefore we need to shift the rhs of the
 * added cut from _alpha_ to _alpha - beta * xbar_, which has taken care of in the routine addCut() in cuts.c. We do not need to worry about the shift
 * of rhs for the dropped cuts.
 * This function performs the change of rhs when the incumbent changes, as described above. */
int changeQPrhs(probType *prob, cellType *cell, vector xk) {
	int 	status = 0, cnt;
	vector 	rhs;
	intvec 	indices;

	if (!(rhs =(vector) arr_alloc(prob->num->rows+cell->cuts->cnt+1, double)))
		errMsg("Allocation", "changeRhs", "rhs",0);
	if (!(indices =(intvec) arr_alloc(prob->num->rows+cell->cuts->cnt, int)))
		errMsg("Allocation", "changeRhs", "indices",0);
	/* Be careful with the one_norm!! In the CxX() routine, it assumes the 0th element is reserved for the 1_norm, in the returned vector, the T sparse
	 vector, and the x vector. */
	for (cnt = 0; cnt < prob->num->rows; cnt++) {
		rhs[cnt + 1] = prob->sp->rhsx[cnt];
		indices[cnt] = cnt;
	}

	/* b - A * xbar */
	rhs = MSparsexvSub(prob->Dbar, xk, rhs);

	/*** new rhs = alpha - beta * xbar (benders cuts)***/
	for (cnt = 0; cnt < cell->cuts->cnt; cnt++) {
		rhs[prob->num->rows+cnt+1] = cell->cuts->vals[cnt]->alpha - vXv(cell->cuts->vals[cnt]->beta, xk, NULL, prob->sp->mac);
		indices[prob->num->rows+cnt] = cell->cuts->vals[cnt]->rowNum;

		cell->cuts->vals[cnt]->alphaIncumb = rhs[prob->num->rows+cnt+1];
	}

	/* Now we change the right-hand of the master problem. */
	status = changeRHS(cell->master->lp, prob->num->rows + cell->cuts->cnt, indices, rhs+1);
	if (status)	{
		errMsg("solver", "changeQPrhs", "failed to change the right-hand side in the solver", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(indices);
	return 0;
}//END changeQPrhs()

/* This function changes the (lower) bounds of the variables, while changing from x to d. The lower bounds of d varibles are -xbar
 * (incumbent solution). */
int changeQPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector xk) {
	int 	status = 0, cnt;
	vector	lbounds, ubounds;
	intvec	lindices, uindices;
	char 	*llu, *ulu;

	if (!(lbounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "changeBounds", "lbounds",0);
	if (!(lindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "change_bounds", "lindices",0);
	if (!(llu = arr_alloc(numCols, char)))
		errMsg("Allocation", "changeBounds", "llu",0);

	if (!(ubounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "change_bounds", "ubounds",0);
	if (!(uindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "changeBounds", "uindices",0);
	if (!(ulu = arr_alloc(numCols, char)))
		errMsg("Allocation", "changeBounds", "ulu",0);

	/* Change the Upper Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		ubounds[cnt] = bdu[cnt] - xk[cnt + 1];
		uindices[cnt] = cnt;
		ulu[cnt] = 'U';
	}

	status = changeBDS(lp, numCols, uindices, ulu, ubounds);
	if (status) {
		errMsg("algorithm", "changeQP", "failed to change the upper bound in the solver", 0);
		return 1;
	}

	/* Change the Lower Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		lbounds[cnt] = bdl[cnt] - xk[cnt + 1];
		lindices[cnt] = cnt;
		llu[cnt] = 'L';
	}

	status = changeBDS(lp, numCols, lindices, llu, lbounds);
	if (status) {
		errMsg("algorithm", "changeQP", "failed to change the lower bound in the solver", 0);
		return 1;
	}

	mem_free(lbounds); mem_free(lindices); mem_free(llu);
	mem_free(ubounds); mem_free(uindices); mem_free(ulu);

	return 0;
}//END changeQPbds()
