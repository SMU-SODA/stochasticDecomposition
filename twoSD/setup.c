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

int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell, vector *meanSol) {
	vector	lb = NULL;
	int 	t;

	/* setup mean value problem which will act as reference for all future computations */
	(*meanSol) = meanProblem(orig, stoc);
	if ( (*meanSol) == NULL ) {
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
	FILE *decompose;
	extern outputDir;
	openFile(outputDir, "decompose.txt", "w");
	printDecomposeSummary(tim, (*prob));
#endif

	/* ensure that we have a linear programs at all stages */
	t = 0;
	while ( t < tim->numStages ) {
		if ( (*prob)[t++]->sp->type  != PROB_LP )
			printf("Warning :: Stage-%d problem is a mixed-integer program. Solving its linear relaxation.\n", t);
	}

	/* create the cells which will be used in the algorithms */
	(*cell) = newCell(stoc, (*prob), (*meanSol));
	if ( (*cell) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to create the necessary cell structure", 0);
		goto TERMINATE;
	}

	mem_free(lb);
	return 0;
	TERMINATE:
	mem_free(lb);
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
	cell->cuts = cell->fCuts = NULL;
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
	cell->gamma 	= 0.0;
	cell->normDk_1 	= 0.0;
	cell->normDk 	= 0.0;

	/* lower bounding approximations held in cuts structure */
	cell->maxCuts = config.CUT_MULT * prob[0]->num->cols + 3;
	cell->cuts 	  = newCuts(cell->maxCuts);

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

	cell->spFeasFlag = TRUE;
	cell->fCuts		= newCuts(cell->maxCuts);
	cell->fCutsPool = newCuts(cell->maxCuts);
	cell->feasCnt 		= 0;
	cell->infeasIncumb 	= FALSE;
	cell->fUpdt[0] = cell->fUpdt[1] = 0;

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
	int 	status, r2 = 1, r3 = 1, maxReps = 30;

	fptr = fopen("config.sd", "r");
	if ( fptr == NULL ) {
		errMsg("read", "readConfig", "failed to open configuration file", 0);
		return 1;
	}
	config.RUN_SEED = (long long *) arr_alloc(maxReps+1, long long);
	config.EVAL_SEED = (long long *) arr_alloc(maxReps+1, long long);
	config.SUBPROB_SAMPLE_SEED = (long long *) arr_alloc(maxReps+1, long long);
	config.NUM_REPS = 0;
	config.MAX_OBS = 0;

	while ((status = (fscanf(fptr, "%s", line) != EOF))) {
		if (!(strcmp(line, "RUN_SEED"))) {
			fscanf(fptr, "%lld", &config.RUN_SEED[config.NUM_REPS+1]);
			config.NUM_REPS++;
			if ( config.NUM_REPS > maxReps ) {
				config.RUN_SEED = (long long *) mem_realloc(config.RUN_SEED, (2*maxReps+1)*sizeof(long long));
				maxReps *= 2;
			}
		}
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
		else if (!(strcmp(line, "EVAL_SEED"))) {
			fscanf(fptr, "%lld", &config.EVAL_SEED[r2++]);
			if ( r2 > maxReps ) {
				config.RUN_SEED = (long long *) mem_realloc(config.RUN_SEED, (2*maxReps+1)*sizeof(long long));
				maxReps *= 2;
			}
		}
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
		else if (!(strcmp(line, "SUBPROB_SAMPLE_SEED"))) {
			fscanf(fptr, "%lld", &config.SUBPROB_SAMPLE_SEED[r3++]);
			if ( r3 > maxReps ) {
				config.RUN_SEED = (long long *) mem_realloc(config.RUN_SEED, (2*maxReps+1)*sizeof(long long));
				maxReps *= 2;
			}
		}
		else if (!(strcmp(line, "SAA")))
			fscanf(fptr, "%d", &config.SAA);

		else if (!(strcmp(line, "MULTIPLE_REP")))
					fscanf(fptr, "%d", &config.MULTIPLE_REP);
		else if (!strcmp(line, "//"))
			fgets(comment, 2*BLOCKSIZE, fptr);
		else {
			printf ("%s\n", line);
			errMsg("read", "readConfig", "unrecognized parameter in configuration file", 1);
		}
	}

	fclose(fptr);
	if ( config.NUM_REPS != (r3 -1) )
		errMsg("read", "readConfig", "number of SUBPROB_SAMPLE_SEED do not match NUM_REPS", 0);

	if ( config.MULTIPLE_REP == 0 )
		config.NUM_REPS = 1;

	return 0;
}//END readConfig()

void freeConfig() {

	if (config.RUN_SEED) mem_free(config.RUN_SEED);
	if (config.EVAL_SEED) mem_free(config.EVAL_SEED);
	if (config.SUBPROB_SAMPLE_SEED) mem_free(config.SUBPROB_SAMPLE_SEED);
}//END freeConfig()

int cleanCellType(cellType *cell, probType *prob, vector xk) {
	int cnt;

	/* constants and arrays */
	cell->k = 0;
	cell->LPcnt = 0;
	cell->optFlag 		 = FALSE;
	cell->dualStableFlag = FALSE;
	cell->spFeasFlag 	 = TRUE;

	copyVector(xk, cell->candidX, prob->num->cols, TRUE);
	cell->candidEst	= prob->lb + vXvSparse(cell->candidX, prob->dBar);

	if (config.MASTER_TYPE == PROB_QP) {
		copyVector(xk, cell->incumbX, prob->num->cols, TRUE);
		cell->incumbEst = cell->candidEst;
		cell->quadScalar= config.MIN_QUAD_SCALAR;
		cell->iCutIdx   = 0;
		cell->iCutUpdt  = 0;
		cell->incumbChg = TRUE;
	}
	cell->gamma 	= 0.0;
	cell->normDk_1 	= 0.0;
	cell->normDk 	= 0.0;

	/* oneProblem structures and solver elements */
	for ( cnt = prob->num->rows+cell->cuts->cnt+cell->fCuts->cnt-1; cnt >= prob->num->rows; cnt-- )
		if (  removeRow(cell->master->lp, cnt, cnt) ) {
			errMsg("solver", "cleanCellType", "failed to remove a row from master problem", 0);
			return 1;
		}
	cell->master->mar = prob->num->rows;
	if( changeQPproximal(cell->master->lp, prob->num->cols, cell->quadScalar)) {
		errMsg("algorithm", "cleanCellType", "failed to change the proximal term", 0);
		return 1;
	}

	/* cuts */
	if (cell->cuts) freeCutsType(cell->cuts, TRUE);
	if (cell->fCuts) freeCutsType(cell->fCuts, TRUE);
	if (cell->fCutsPool) freeCutsType(cell->fCutsPool, TRUE);
	cell->feasCnt 		= 0;
	cell->infeasIncumb 	= FALSE;
	cell->fUpdt[0] = cell->fUpdt[1] = 0;

	/* stochastic components */
	if (cell->basis) freeBasisType(cell->basis, TRUE);
	if (cell->delta) freeDeltaType(cell->delta, cell->lambda->cnt, cell->omega->cnt, TRUE);
	if (cell->lambda) freeLambdaType(cell->lambda, TRUE);
	if (cell->sigma) freeSigmaType(cell->sigma, TRUE);
	if (cell->omega) freeOmegaType(cell->omega, TRUE);

	if ( config.MASTER_TYPE == PROB_QP ) {
		if ( constructQP(prob, cell, cell->incumbX, cell->quadScalar) ) {
			errMsg("setup", "newCell", "failed to change the right-hand side after incumbent change", 0);
			return NULL;
		}

		cell->incumbChg = FALSE;
#if defined(SETUP_CHECK)
		if ( writeProblem(cell->aster->lp, "cleanedQPMaster.lp") ) {
			errMsg("write problem", "new_master", "failed to write master problem to file",0);
			return NULL;
		}
#endif
	}

	return 0;
}//END cleanCellType()

void freeCellType(cellType *cell) {

	if ( cell ) {
		if (cell->master) freeOneProblem(cell->master);
		if (cell->candidX) mem_free(cell->candidX);
		if (cell->incumbX) mem_free(cell->incumbX);
		if (cell->piM) mem_free(cell->piM);
		if (cell->djM) mem_free(cell->djM);
		if (cell->cuts) freeCutsType(cell->cuts, FALSE);
		if (cell->fCuts) freeCutsType(cell->fCuts, FALSE);
		if (cell->fCutsPool) freeCutsType(cell->fCutsPool, FALSE);
		if (cell->basis) freeBasisType(cell->basis, FALSE);
		if (cell->delta) freeDeltaType(cell->delta, cell->lambda->cnt, cell->omega->cnt, FALSE);
		if (cell->omega) freeOmegaType(cell->omega, FALSE);
		if (cell->lambda) freeLambdaType(cell->lambda, FALSE);
		if (cell->sigma) freeSigmaType(cell->sigma, FALSE);
		if (cell->pi_ratio) mem_free(cell->pi_ratio);
		mem_free(cell);
	}

}//END freeCellType()
