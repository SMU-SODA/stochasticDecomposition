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
	config.SAMPLE_INCREMENT = 1;

	while ((status = (fscanf(fptr, "%s", line) != EOF))) {
		if (!(strcmp(line, "RUN_SEED"))) {
			fscanf(fptr, "%lld", &config.RUN_SEED[config.NUM_SEEDS+1]);
			config.NUM_SEEDS++;
			if ( config.NUM_SEEDS > maxReps ) {
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
		else if (!(strcmp(line, "MAX_ITER_CLBK")))
			fscanf(fptr, "%d", &config.MAX_ITER_CLBK);
		else if (!(strcmp(line, "MAX_ITER_ROOT")))
			fscanf(fptr, "%d", &config.MAX_ITER_ROOT);
		else if (!(strcmp(line, "MAX_NODES")))
			fscanf(fptr, "%d", &config.MAX_NODES);
		else if (!(strcmp(line, "MASTER_TYPE")))
			fscanf(fptr, "%d", &config.MASTER_TYPE);
		else if (!(strcmp(line, "VAR_STR")))
			fscanf(fptr, "%d", &config.VAR_STR);
		else if (!(strcmp(line, "HEURST_FLAG")))
			fscanf(fptr, "%d", &config.HEURST_FLAG);
		else if (!(strcmp(line, "BRN_STR")))
			fscanf(fptr, "%d", &config.BRN_STR);
		else if (!(strcmp(line, "Pi_EVAL_FLAG")))
			fscanf(fptr, "%d", &config.Pi_EVAL_FLAG);
		else if (!(strcmp(line, "CUT_MULT")))
			fscanf(fptr, "%d", &config.CUT_MULT);
		else if (!(strcmp(line, "NodeNum")))
			fscanf(fptr, "%d", &config.NodeNum);
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
		else if (!(strcmp(line, "SMIP_SOLVER")))
			fscanf(fptr, "%d", &config.SMIP);
		else if (!(strcmp(line, "SMIP_ALGO")))
			fscanf(fptr, "%d", &config.ALGO);
		else if (!(strcmp(line, "SMIP_OPTGAP")))
			fscanf(fptr, "%lf", &config.SMIP_OPTGAP);
		else if (!(strcmp(line, "DUAL_STABILITY")))
			fscanf(fptr, "%d", &config.DUAL_STABILITY);
		else if (!(strcmp(line, "PI_EVAL_START")))
			fscanf(fptr, "%d", &config.PI_EVAL_START);
		else if (!(strcmp(line, "PI_CYCLE")))
			fscanf(fptr, "%d", &config.PI_CYCLE);
		else if (!(strcmp(line, "PERCENT_PASS")))
			fscanf(fptr, "%lf", &config.PERCENT_PASS);
		else if (!(strcmp(line, "MIN_X")))
			fscanf(fptr, "%lf", &config.MIN_X);
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
		else if (!(strcmp(line, "INT_TOLERANCE")))
			fscanf(fptr, "%lf", &config.INT_TOLERANCE);
		else if (!strcmp(line, "//"))
			fgets(comment, 2*BLOCKSIZE, fptr);
		else {
			printf ("%s\n", line);
			errMsg("read", "readConfig", "unrecognized parameter in configuration file", 1);
		}
	}

	fclose(fptr);

	config.NUM_SEEDS = minimum(config.NUM_SEEDS, r2);
	if ( config.MULTIPLE_REP > config.NUM_SEEDS ) {
		printf("Requesting to perform more replications than the number of seeds provided.\n");
		return 1;
	}

	if ( config.MULTIPLE_REP == 1 ) {
		config.COMPROMISE_PROB = 0;
	}

	return 0;
}//END readConfig()

int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell,
		batchSummary **batch, dVector *meanSol, int type) {
	int 	t;
	dVector lb = NULL;

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
		mem_free(lb);
		return 1;
	}


	/* decompose the problem into master and subproblem */
	(*prob) = newProb(orig, stoc, tim, lb, config.TOLERANCE);
	if ((*prob) == NULL) {
		errMsg("setup", "setupAlgo", "failed to update probType with elements specific to algorithm", 0);
		mem_free(lb);
		return 1;
	}

#ifdef DECOMPOSE_CHECK
	printDecomposeSummary(stdout, orig->name, tim, (*prob));
#endif

	/* ensure that we have a linear programs at all stages */
	t = 0;
	while (t < tim->numStages) {
		if (config.SMIP != 0)
		{
			if ((*prob)[t++]->sp->type != PROB_LP)
				printf("Warning :: Stage-%d problem is a mixed-integer program. Solving its linear relaxation.\n", t);
		}
		else
		{
			t++;
		}
	}


	/* create the cells which will be used in the algorithms */
	(*cell) = newCell(stoc, (*prob), (*meanSol), type);
	if ( (*cell) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to create the necessary cell structure", 0);
		mem_free(lb);
		return 1;
	}

	if ( config.MULTIPLE_REP > 1 )
		(*batch)  = newBatchSummary((*prob)[0], config.MULTIPLE_REP);

	mem_free(lb);
	return 0;
}//END setupAlgo()

 /* creating info summary from bnc */
bncInfoSummary *setupBncInfo(int cnt) {
	bncInfoSummary *info = NULL;

	if (!(info = (bncInfoSummary *) mem_malloc(sizeof(bncInfoSummary))))
		errMsg("Memory allocation", "newCell", "failed to allocate memory to bncInfo", 0);
	if (!(info->disjnct = (iVector)arr_alloc(cnt, int)))
		errMsg("allocation", "setupBncInfo", "info->disjnct", 0);
	if (!(info->vals = (dVector)arr_alloc(cnt, double)))
		errMsg("allocation", "setupBncInfo", "info->vals", 0);

	return info;
}//END setupBncInfo()

/* This function is used to create cells used in the algorithm */
cellType *newCell(stocType *stoc, probType **prob, dVector xk, int type) {
	cellType    *cell = NULL;
	int			length;

	/* allocate memory to all cells used in the algorithm. The first cell belongs to the master problem, while the rest correspond to each of the
	 * sub-agents in the problem.  */
	cell = (cellType *) mem_malloc(sizeof(cellType));

	cell->master = cell->subprob = NULL;
	cell->candidX = cell->incumbX = NULL;
	cell->piM = cell->djM = NULL;
	cell->activeCuts = cell->fcuts = NULL;
	cell->cutsPool = NULL;
	cell->lambda = NULL; cell->sigma = NULL; cell->delta = NULL; cell->omega = NULL;
	cell->pi_ratio = NULL;

	/* setup the master problem */
	cell->master = newMaster(prob[0]->sp, prob[0]->lb, type);
	if (cell->master == NULL) {
		errMsg("setup", "newCell", "failed to setup the master problem", 0);
		return NULL;
	}
	
	/* setup the subproblem */
	cell->subprob = newSubprob(prob[1]->sp);
	
	/* setup the BnB part of the cell */
	cell->cur_rowname = NULL;
	cell->isinBnB = false;

	/* -+-+-+-+-+-+-+-+-+-+-+ Allocating memory to other variables that belongs to master mcell +-+-+-+-+-+-+-+-+-+- */
	cell->k 	= 0;
	cell->ki    = 0;
	cell->mk    = 0;
	cell->gk    = 0;
	cell->sampleSize = 0;
	cell->LPcnt = 0;
	if (prob[0]->lb == 0)
		cell->lbType = TRIVIAL;
	else
		cell->lbType = NONTRIVIAL;
	cell->lb = prob[0]->lb;

	/* candidate solution and estimates */
	cell->candidX 			= duplicVector(xk, prob[0]->num->cols, true);
	cell->candidEst 		= prob[0]->lb + vXvSparse(cell->candidX, prob[0]->dBar);

	/* incumbent solution and estimates */
	if ( cell->master->type == PROB_QP ) {
		cell->incumbX   = duplicVector(xk, prob[0]->num->cols, true);
		cell->incumbEst = cell->candidEst;
		cell->quadScalar= config.MIN_QUAD_SCALAR;     						/* The quadratic scalar, 'sigma'*/
		cell->iCutIdx   = 0;
		cell->iCutUpdt  = 0;
		cell->incumbChg = true;
	}
	else {
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
	cell->maxMIPCuts = 1000;
	cell->etaIdx = prob[0]->num->cols;
	cell->activeCuts = newCuts(cell->maxCuts);

	cell->cutsPool = (cutsType **) arr_alloc(config.NodeNum, cutsType *);
	cell->cutsPool[0] = newCuts(cell->maxCuts);
	cell->numPools = 1;

	cell->MIRcuts = newCuts(cell->maxMIPCuts);
	cell->GMIcuts = newCuts(cell->maxMIPCuts);

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

	cell->optFlag = false;
	cell->MIPFlag = false;

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
	cell->masterFeasFlag = true;
	cell->fcuts		= newCuts(cell->maxCuts);
	cell->fcutsPool = newCuts(cell->maxCuts);
	cell->feasCnt 		= 0;
	cell->infeasIncumb 	= false;
	cell->fUpdt[0] = cell->fUpdt[1] = 0;

	cell->time.repTime = cell->time.iterTime = cell->time.masterIter = cell->time.subprobIter = cell->time.optTestIter = cell->time.argmaxIter = 0.0;
	cell->time.iterAccumTime = cell->time.masterAccumTime = cell->time.subprobAccumTime = cell->time.optTestAccumTime = cell->time.argmaxAccumTime = 0.0;

	/* construct the QP using the current incumbent */
	if ( cell->master->type == PROB_QP ) {
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

	cell->tot_nodes = 0;
	cell->depth = 0;
	cell->d_nodes = 0;
	cell->int_nodes = 0;
	cell->maxiter_nodes = 0;

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
	cell->ki = 0;
	cell->numPools = 0;
	cell->sampleSize = 0;
	cell->LPcnt = 0;
	cell->optFlag 		 = false;
	cell->spFeasFlag 	 = true;
	if ( config.DUAL_STABILITY )
		cell->dualStableFlag 	= false;

	copyVector(xk, cell->candidX, prob->num->cols, true);
	cell->candidEst	= prob->lb + vXvSparse(cell->candidX, prob->dBar);

	if ( cell->master->type == PROB_QP ) {
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

	/* Get the total number of rows */
	int row_num = getNumRows(cell->master->lp);

	/* oneProblem structures and solver elements */
	if (prob->num->rows < row_num)
	{
		for (cnt = row_num - 1; cnt >= prob->num->rows; cnt--)
			if (removeRow(cell->master->lp, cnt, cnt)) {
				printf("row Num %d - tot rows %d - orig rows %d", cnt, row_num, prob->num->rows);
				errMsg("solver", "cleanCellType", "failed to remove a row from master problem", 0);
				return 1;
			}
	}

	cell->master->mar = prob->num->rows;
	if( changeQPproximal(cell->master->lp, prob->num->cols, cell->quadScalar)) {
		errMsg("algorithm", "cleanCellType", "failed to change the proximal term", 0);
		return 1;
	}

	/* cuts */
	if (cell->activeCuts) freeCutsType(cell->activeCuts, true);
	if (cell->MIRcuts) freeCutsType(cell->MIRcuts, true);
	if (cell->GMIcuts) freeCutsType(cell->GMIcuts, true);
	if (cell->fcuts) freeCutsType(cell->fcuts, true);
	if (cell->fcutsPool) freeCutsType(cell->fcutsPool, true);
	if (*cell->cutsPool) {
		for (int n = 0; n < (*cell->cutsPool)->cnt; n++)
			if (cell->cutsPool[n]) freeCutsType(cell->cutsPool[n], true);
	}
	cell->cutsPool[0] = newCuts(cell->maxCuts);
	cell->numPools = 1;
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

	if ( cell->master->type == PROB_QP ) {
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
		if (cell->activeCuts) freeCutsType(cell->activeCuts, false);
		if (cell->MIRcuts) freeCutsType(cell->MIRcuts, false);
		if (cell->GMIcuts) freeCutsType(cell->GMIcuts, false);
		if (cell->fcuts) freeCutsType(cell->fcuts, false);
		if (cell->fcutsPool) freeCutsType(cell->fcutsPool, false);
		if (cell->delta) freeDeltaType(cell->delta, cell->lambda->cnt, cell->omega->cnt, false);
		if (cell->omega) freeOmegaType(cell->omega, false);
		if (cell->lambda) freeLambdaType(cell->lambda, false);
		if (cell->sigma) freeSigmaType(cell->sigma, false);
		if (cell->basis) freeBasisType(cell->basis, false);
		if (cell->sample) freeSampleType(cell->sample);
		if (cell->pi_ratio) mem_free(cell->pi_ratio);
		if (*cell->cutsPool) {
			for ( int n = 0; n < (*cell->cutsPool)->cnt; n++ )
				if ( cell->cutsPool[n] ) freeCutsType(cell->cutsPool[n], false);
			mem_free(cell->cutsPool);
		}
		mem_free(cell);
	}

}//END freeCellType()