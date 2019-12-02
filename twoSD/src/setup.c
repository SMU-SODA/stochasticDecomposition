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

#include "twoSD.h"
#include <stdio.h>

extern configType config;

int readConfig(cString path2config, cString inputDir) {
	FILE 	*fptr;
	char	line[2*BLOCKSIZE], comment[2*BLOCKSIZE];
	int 	status, r2 = 1, maxReps = 30;

	strcpy(line, path2config);	strcat(line, "config.sd");
	fptr = fopen(line, "r");
	if ( fptr == NULL ) {
		errMsg("read", "readConfig", "failed to open configuration file", 0);
		return 1;
	}

	config.RUN_SEED = (long long *) arr_alloc(maxReps+1, long long);
	config.EVAL_SEED = (long long *) arr_alloc(maxReps+1, long long);
	config.NUM_REPS = 0;
	config.SAMPLE_INCREMENT = 1;

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
		else if (!(strcmp(line, "CUT_MULT")))
			fscanf(fptr, "%d", &config.CUT_MULT);
		else if (!(strcmp(line, "TAU")))
			fscanf(fptr, "%d", &config.TAU);
		else if (!(strcmp(line, "MIN_QUAD_SCALAR")))
			fscanf(fptr, "%lf", &config.MIN_QUAD_SCALAR);
		else if (!(strcmp(line, "MAX_QUAD_SCALAR")))
			fscanf(fptr, "%lf", &config.MAX_QUAD_SCALAR);
		else if (!(strcmp(line, "R1")))
			fscanf(fptr, "%lf", &config.R1);
		else if (!(strcmp(line, "R2")))
			fscanf(fptr, "%lf", &config.R2);
		else if (!(strcmp(line, "R3")))
			fscanf(fptr, "%lf", &config.R3);
		else if (!(strcmp(line, "DUAL_STABILITY")))
			fscanf(fptr, "%d", &config.DUAL_STABILITY);
		else if (!(strcmp(line, "PI_EVAL_START")))
			fscanf(fptr, "%d", &config.PI_EVAL_START);
		else if (!(strcmp(line, "PI_CYCLE")))
			fscanf(fptr, "%d", &config.PI_CYCLE);
		else if (!(strcmp(line, "PERCENT_PASS")))
			fscanf(fptr, "%lf", &config.PERCENT_PASS);
		else if (!(strcmp(line, "SCAN_LEN")))
			fscanf(fptr, "%d", &config.SCAN_LEN);
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
		else if (!(strcmp(line, "PRE_EPSILON")))
			fscanf(fptr, "%lf", &config.PRE_EPSILON);
		else if (!(strcmp(line, "EPSILON")))
			fscanf(fptr, "%lf", &config.EPSILON);
		else if (!(strcmp(line, "BOOTSTRAP_REP")))
			fscanf(fptr, "%d", &config.BOOTSTRAP_REP);

		else if (!(strcmp(line, "MULTIPLE_REP")))
			fscanf(fptr, "%d", &config.MULTIPLE_REP);
		else if (!(strcmp(line, "COMPROMISE_PROB")))
			fscanf(fptr, "%d", &config.COMPROMISE_PROB);

		else if (!(strcmp(line, "SAMPLE_INCREMENT")))
			fscanf(fptr, "%d", &config.SAMPLE_INCREMENT);

		else if (!strcmp(line, "//"))
			fgets(comment, 2*BLOCKSIZE, fptr);
		else {
			printf ("%s\n", line);
			errMsg("read", "readConfig", "unrecognized parameter in configuration file", 1);
		}
	}

	fclose(fptr);

	if ( config.MULTIPLE_REP == 0 ) {
		config.NUM_REPS = 1;
		config.COMPROMISE_PROB = 0;
	}

	return 0;
}//END readConfig()

int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell,
		batchSummary **batch, dVector *meanSol) {
	dVector	lb = NULL;
	int 	t;

	/* setup mean value problem which will act as reference for all future computations */
	(*meanSol) = meanProblem(orig, stoc);
	if ( meanSol == NULL ) {
		errMsg("setup", "setupAlgo", "failed to setup and solve mean value problem", 0);
		return 1;
	}

	/* calculate lower bounds for each stage */
	lb = calcLowerBound(orig, tim, stoc);
	if ( lb == NULL )  {
		errMsg("setup", "setupAlgo", "failed to compute lower bounds on stage problem", 0);
		mem_free(lb); return 1;
	}

	/* decompose the problem into master and subproblem */
	(*prob) = newProb(orig, stoc, tim, lb, config.TOLERANCE);
	if ( (*prob) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to update probType with elements specific to algorithm", 0);
		mem_free(lb); return 1;
	}

#ifdef DECOMPOSE_CHECK
	printDecomposeSummary(stdout, orig->name, tim, (*prob));
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
		return 1;
	}

	if ( config.NUM_REPS > 1 )
		(*batch)  = newBatchSummary((*prob)[0], config.NUM_REPS);

	mem_free(lb);
	return 0;
}//END setupAlgo()

/* This function is used to create cells used in the algorithm */
cellType *newCell(stocType *stoc, probType **prob, dVector xk) {
	cellType    *cell = NULL;
	int			length;

	/* allocate memory to all cells used in the algorithm. The first cell belongs to the master problem, while the rest correspond to each of the
	 * sub-agents in the problem.  */
	if (!(cell = (cellType *) mem_malloc(sizeof(cellType))) )
		errMsg("Memory allocation", "newCell", "failed to allocate memory to cell",0);
	cell->master = cell->subprob = NULL;
	cell->candidX = cell->incumbX = NULL;
	cell->piM = cell->djM = NULL;
	cell->cuts = cell->fcuts = NULL;
	cell->lambda = NULL; cell->sigma = NULL; cell->delta = NULL; cell->omega = NULL;
	cell->pi_ratio = NULL;

	/* setup the master problem */
	cell->master = newMaster(prob[0]->sp, prob[0]->lb);
	if ( cell->master == NULL ) {
		errMsg("setup", "newCell", "failed to setup the master problem", 0);
		return NULL;
	}
	/* setup the subproblem */
	cell->subprob = newSubprob(prob[1]->sp);

	/* -+-+-+-+-+-+-+-+-+-+-+ Allocating memory to other variables that belongs to master mcell +-+-+-+-+-+-+-+-+-+- */
	cell->k 	= 0;
	cell->sampleSize = 0;
	cell->LPcnt = 0;
	if (prob[0]->lb == 0)
		cell->lbType = TRIVIAL;
	else
		cell->lbType = NONTRIVIAL;
	cell->lb = prob[0]->lb;

	/* candidate solution and estimates */
	cell->candidX 			= duplicVector(xk, prob[0]->num->cols);
	cell->candidEst 		= prob[0]->lb + vXvSparse(cell->candidX, prob[0]->dBar);

	/* incumbent solution and estimates */
	if (config.MASTER_TYPE == PROB_QP) {
		cell->incumbX   = duplicVector(xk, prob[0]->num->cols);
		cell->incumbMIPX = duplicVector(xk, prob[0]->num->cols);
		cell->incumbEst = cell->candidEst;
		cell->quadScalar= config.MIN_QUAD_SCALAR;     						/* The quadratic scalar, 'sigma'*/
		cell->iCutIdx   = 0;
		cell->iCutUpdt  = 0;
		cell->incumbChg = true;
	}
	else {
		cell->incumbMIPX = NULL;
		cell->incumbX   = NULL;
		cell->incumbEst = 0.0;
		cell->quadScalar= 0.0;
		cell->iCutIdx   = -1;
		cell->iCutUpdt  = -1;
		cell->incumbChg = false;
	}
	cell->gamma 			= 0.0;
	cell->normDk_1 			= 0.0;
	cell->normDk 			= 0.0;

	/* lower bounding approximations held in cuts structure */
	cell->maxCuts = config.CUT_MULT * prob[0]->num->cols + 3;
	cell->cuts 	  = newCuts(cell->maxCuts);
	cell->MIPcuts = newCuts(cell->maxCuts);

	/* solution parts of the cell */
	if ( !(cell->djM = (dVector) arr_alloc(prob[0]->num->cols + 2, double)) )
		errMsg("allocation", "newCell", "cell->di", 0);
	if ( !(cell->piM = (dVector) arr_alloc(prob[0]->num->rows + cell->maxCuts + 1, double)) )
		errMsg("allocation", "newCell", "cell->piM", 0);

	/* stochastic elements: we need more room to store basis information when the cost coefficients are random. */
	if ( prob[1]->num->rvdOmCnt > 0 )
		length = prob[1]->num->rvdOmCnt*config.MAX_ITER + config.MAX_ITER / config.TAU + 1;
	else
		length = config.MAX_ITER + config.MAX_ITER / config.TAU + 1;
	cell->basis  = newBasisType(config.MAX_ITER, prob[1]->num->cols, prob[1]->num->rows, WORDLENGTH);
	cell->lambda = newLambda(length, 0, prob[1]->num->rvRowCnt);
	cell->sigma  = newSigma(length, prob[1]->num->cntCcols, 0);
	cell->delta  = newDelta(length);
	cell->omega  = newOmega(prob[1]->num->numRV, config.MAX_ITER);
	cell->sample = newSample(config.SAMPLE_INCREMENT);

	cell->optFlag 			= false;

	/* Dual stability test is disabled DUAL_STABILITY is false. */
	if ( !config.DUAL_STABILITY ) {
		cell->dualStableFlag = true;
		cell->pi_ratio = NULL;
	}
	else {
		cell->dualStableFlag 	= false;
		if ( !(cell->pi_ratio = (dVector) arr_alloc(config.SCAN_LEN, double)) )
			errMsg("allocation", "newCell", "cell->pi_ratio", 0);
	}

	cell->spFeasFlag = true;
	cell->fcuts		= newCuts(cell->maxCuts);
	cell->fcutsPool = newCuts(cell->maxCuts);
	cell->feasCnt 		= 0;
	cell->infeasIncumb 	= false;
	cell->fUpdt[0] = cell->fUpdt[1] = 0;

	cell->time.repTime = cell->time.iterTime = cell->time.masterIter = cell->time.subprobIter = cell->time.optTestIter = cell->time.argmaxIter = 0.0;
	cell->time.iterAccumTime = cell->time.masterAccumTime = cell->time.subprobAccumTime = cell->time.optTestAccumTime = cell->time.argmaxAccumTime = 0.0;

	/* construct the QP using the current incumbent */
	if ( config.MASTER_TYPE == PROB_QP ) {
		if ( constructQP(prob[0], cell, cell->incumbX, cell->quadScalar) ) {
			errMsg("setup", "newCell", "failed to change the right-hand side after incumbent change", 0);
			return NULL;
		}

		cell->incumbChg = false;
#if defined(SETUP_CHECK)
		if ( writeProblem(cell->aster->lp, "newQPMaster.lp") ) {
			errMsg("write problem", "new_master", "failed to write master problem to file",0);
			return NULL;
		}
#endif
	}

	return cell;
}//END newCell()

void freeConfig() {

	if (config.RUN_SEED) mem_free(config.RUN_SEED);
	if (config.EVAL_SEED) mem_free(config.EVAL_SEED);

}//END freeConfig()

int cleanCellType(cellType *cell, probType *prob, dVector xk) {
	int cnt;

	/* constants and arrays */
	cell->k = 0;
	cell->sampleSize = 0;
	cell->LPcnt = 0;
	cell->optFlag 		 = false;
	cell->spFeasFlag 	 = true;
	if ( config.DUAL_STABILITY )
		cell->dualStableFlag 	= false;

	copyVector(xk, cell->candidX, prob->num->cols, true);
	cell->candidEst	= prob->lb + vXvSparse(cell->candidX, prob->dBar);

	if (config.MASTER_TYPE == PROB_QP) {
		copyVector(xk, cell->incumbX, prob->num->cols, true);
		cell->incumbEst = cell->candidEst;
		cell->quadScalar= config.MIN_QUAD_SCALAR;
		cell->iCutIdx   = 0;
		cell->iCutUpdt  = 0;
		cell->incumbChg = true;
	}
	cell->gamma 	= 0.0;
	cell->normDk_1 	= 0.0;
	cell->normDk 	= 0.0;

	/* oneProblem structures and solver elements */
	for ( cnt = prob->num->rows+cell->cuts->cnt+cell->fcuts->cnt-1; cnt >= prob->num->rows; cnt-- )
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
	if (cell->cuts) freeCutsType(cell->cuts, true);
	if (cell->fcuts) freeCutsType(cell->fcuts, true);
	if (cell->fcutsPool) freeCutsType(cell->fcutsPool, true);
	cell->feasCnt 		= 0;
	cell->infeasIncumb 	= false;
	cell->fUpdt[0] = cell->fUpdt[1] = 0;

	/* stochastic components */
	if (cell->basis) freeBasisType(cell->basis, true);
	if (cell->delta) freeDeltaType(cell->delta, cell->lambda->cnt, cell->omega->cnt, true);
	if (cell->lambda) freeLambdaType(cell->lambda, true);
	if (cell->sigma) freeSigmaType(cell->sigma, true);
	if (cell->omega) freeOmegaType(cell->omega, true);

	/* reset all the clocks */
	cell->time.repTime = cell->time.iterTime = cell->time.masterIter = cell->time.subprobIter = cell->time.optTestIter = cell->time.argmaxIter = 0.0;
	cell->time.iterAccumTime = cell->time.masterAccumTime = cell->time.subprobAccumTime = cell->time.optTestAccumTime = cell->time.argmaxAccumTime = 0.0;

	if ( config.MASTER_TYPE == PROB_QP ) {
		if ( constructQP(prob, cell, cell->incumbX, cell->quadScalar) ) {
			errMsg("setup", "newCell", "failed to change the right-hand side after incumbent change", 0);
			return 1;
		}

		cell->incumbChg = false;
#if defined(SETUP_CHECK)
		if ( writeProblem(cell->aster->lp, "cleanedQPMaster.lp") ) {
			errMsg("write problem", "new_master", "failed to write master problem to file",0);
			return 1;
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
		if (cell->cuts) freeCutsType(cell->cuts, false);
		if (cell->fcuts) freeCutsType(cell->fcuts, false);
		if (cell->fcutsPool) freeCutsType(cell->fcutsPool, false);
		if (cell->delta) freeDeltaType(cell->delta, cell->lambda->cnt, cell->omega->cnt, false);
		if (cell->omega) freeOmegaType(cell->omega, false);
		if (cell->lambda) freeLambdaType(cell->lambda, false);
		if (cell->sigma) freeSigmaType(cell->sigma, false);
		if (cell->basis) freeBasisType(cell->basis, false);
		if (cell->sample) freeSampleType(cell->sample);
		if (cell->pi_ratio) mem_free(cell->pi_ratio);
		mem_free(cell);
	}

}//END freeCellType()
