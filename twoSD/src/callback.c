/*
 * callback.c
 *
 *  Created on: Apr 3, 2020
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"

extern cString outputDir;

//ENVptr	env;
typedef struct {
	int	call;
	int MIPcuts;
	int OPTcuts;
	probType **prob;
	cellType *cell;
	cellType *orig_cell;
	stocType *stoc;
} callbackArgs;

int bendersCallback(stocType *stoc, probType **prob, cellType *cell);
static int CPXPUBLIC usersolve (CPXCENVptr env, void *cbdata, int wherefrom, callbackArgs *args);

extern configType config;

int bendersCallback(stocType *stoc, probType **prob, cellType *cell) {
	iVector indices;
	callbackArgs * cbhandle;
	int status = 0;

	printf("\nStarting branch-and-cut phase of the algorithm....\n");
	printf("Iteration-%4d: ", cell->k); fflush(stdout);
	cell->callback = true;


#if defined(CALLBACK_CHECK)
	writeProblem(cell->master->lp, "callback_init.lp");
#endif

	/* Set up to use MIP callback structure and function. */
	cbhandle = (callbackArgs *)mem_malloc(sizeof(callbackArgs));
	cbhandle->cell = cell; cbhandle->prob = prob; cbhandle->call = 0; cbhandle->stoc = stoc;
	cbhandle->MIPcuts = 0; cbhandle->OPTcuts = getNumRows(cell->master->lp);

	/* Setup the callback solver function */
	if (setsolvecallbackfunc(usersolve, cbhandle)) { return 1; }

	/* update the maximum number of SD iteration */
	config.MAX_ITER += config.MAX_ITER;


	setIntParam(CPXPARAM_Preprocessing_Linear, CPX_OFF);		/* Assure linear mappings between the presolved and original models */
	setIntParam(CPXPARAM_Preprocessing_Presolve, CPX_OFF);		/* Assure linear mappings between the presolved and original models */
	setIntParam(CPX_PARAM_DEPIND, CPX_OFF);				     	/* redundancy checker */
	setIntParam(CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);			/* Let MIP callbacks work on the original model */
	setIntParam(CPXPARAM_MIP_Interval, 1);									/* Set MIP log interval to 1 */
	setIntParam(CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL);	/* Turn on traditional search for use with control callbacks */
	setIntParam(CPX_PARAM_THREADS, 1);										/* Set the number of threads to be used to one */
	setIntParam(CPX_PARAM_HEURFREQ, -1);										/* Turn-off the heuristics */
	setIntParam(CPX_PARAM_GUBCOVERS, -1);								    /* Turn-off the generalized upper bound (GUB) cover cuts */
	setIntParam(CPX_PARAM_FLOWCOVERS, -1);								    /* Turn-off flow cover cuts */
	setIntParam(CPX_PARAM_FLOWPATHS, -1);								    /* Turn-off flow path cuts */
	setIntParam(CPX_PARAM_IMPLBD, -1);									    /* Turn-off globally valid implied bound cuts */
	setIntParam(CPX_PARAM_LANDPCUTS, -1);								    /* Turn-off lift-and-project cuts */
	setIntParam(CPX_PARAM_MCFCUTS, -1);								        /* Turn-off multi-commodity flow cuts */
#if !defined(UserMIPcutsActive)
#if defined(CpxMIRCutsActive)
	setIntParam(CPX_PARAM_MIRCUTS, 1);								        /* Turn-off MIR cuts (mixed integer rounding cuts) */
#else
	setIntParam(CPX_PARAM_MIRCUTS, -1);								        /* Turn-off MIR cuts (mixed integer rounding cuts) */
#endif // defined(CpxMIRCutsActive)

#if defined(CpxGMICutsActive)
	setIntParam(CPX_PARAM_FRACCUTS, 1);								        /* Turn-off Gomory fractional cuts */
#else
	setIntParam(CPX_PARAM_FRACCUTS, -1);								    /* Turn-off Gomory fractional cuts */
#endif // defined(CpxGMICutsActive)
#else
	setIntParam(CPX_PARAM_FRACCUTS, -1);								    /* Turn-off Gomory fractional cuts */
	setIntParam(CPX_PARAM_MIRCUTS, -1);								        /* Turn-off MIR cuts (mixed integer rounding cuts) */
#endif // !defined(UserMIPcutsActive)

#if defined(CALLBACK_WRITE_LP)
	char fname[NAMESIZE];
	sprintf(fname, "%s_%d.lp", "callback_before", 0);
	writeProblem(cell->master->lp, fname);
#endif // defined(CALLBACK_WRITE_LP)

	cell->nodeSol[0].sol = duplicVector(cell->candidX, prob[0]->num->cols);
	cell->nodeSol[0].sol_size = prob[0]->num->cols;
	cell->nodeSol[0].nodeNum = 0;
	cell->nodeSol[0].LB = cell->candidEst;
	cell->nodeSol[0].piM = duplicVector(cell->piM, cell->master->mar);
	cell->nodeSol[0].djM = duplicVector(cell->djM, cell->master->mar);
	cell->nodeSol[0].mar = cell->master->mar;
	cell->nodeSol[0].isInt = isInteger(cell->candidX, prob[0]->num->cols,
		                               0, prob[0]->num->cols - 1, config.INT_TOLERANCE);

	/* Launch the solver in callback mode to solve the master problem */
	if (solveProblem(cell->master->lp, cell->master->name, cell->master->type, cell->master->mar, cell->master->mac,
		&status, config.SMIP_OPTGAP)) {
		errMsg("algorithm", "bendersCallback", "failed to solve the master problem", 0);
		return 1;
	}

	/*Print the summary of callbacks*/
	printf("\ncallback summary:\n");
	printf("\# of callbacks calls:%i  - # of OPT cuts:%i - # of MIP cuts:%i\n", cbhandle->call, cbhandle->OPTcuts, cbhandle->MIPcuts);
	printNodeInfo(cell->nodeSol, cbhandle->call);
	copyVector(cell->candidX, cell->incumbX, prob[0]->num->cols, 1);
	cell->callback = false;
	mem_free(cbhandle);
	return 0;
}//END bendersCallback()

static int CPXPUBLIC usersolve (CPXCENVptr env, void *cbdata, int wherefrom, callbackArgs *args) {
	iVector	candidCuts;
	bool breakLoop = false;
	dVector observ;
	clock_t	tic;
	int status;

//#if defined(CALLBACK_WRITE_LP)
//	char fname[80];
//	sprintf(fname, "%s_%d.lp", "clbk_beforelp", args->call);
//	writeProblem(args->cell->master->lp, fname);
//#endif // defined(CALLBACK_WRITE_LP)

	observ = (dVector)arr_alloc(args->stoc->numOmega + 1, double);
	args->cell->ki++;
	args->cell->kii = 1;

	//if ( config.MULTICUT)
	if (false)
		candidCuts = (iVector) arr_alloc(args->cell->omega->cnt, int);
	else
		candidCuts = (iVector) arr_alloc(1, int);

	//int numCol = getNumCols(args->cell->master->lp);
	/* Save the LP pointer from the original master problem. */
	LPptr nodelp = args->cell->master->lp; LPptr temp;
	args->cell->master->lp = NULL;

	/* Get pointer to LP subproblem */
	if ( getcallbacknodelp (cbdata, wherefrom, &temp) ) { return 1; }
	args->cell->master->lp = temp;

#if defined(CALLBACK_WRITE_LP)
	char fname[80];
	sprintf(fname, "%s_%d.lp", "clbk_afterlp", args->call);
	writeProblem(args->cell->master->lp, fname);
#endif // defined(CALLBACK_WRITE_LP)

	/* Get the number of column because if the a variable is deleted
	   the index of \eta variable should be updated */
	int numCol = getNumCols(temp);

	/*Update couns of callback calls*/
	args->call++;

#if defined(CALLBACK_CHECK)
	char fname[NAMESIZE];
	printf("\tCallback routines execution #: %d\n", args->call++);
	sprintf(fname, "%s_%d.lp", "callback", args->call);
	writeProblem(args->cell->master->lp, fname);
	callbackNodeSummary(cbdata, wherefrom, args->cell->master->lp);
#endif


	if ( solveProblem(args->cell->master->lp, args->cell->master->name, PROB_LP, args->cell->master->mar, args->cell->master->mac,
			&status, config.SMIP_OPTGAP) ) {
		writeProblem(args->cell->master->lp, "error.lp");
		errMsg("algorithm", "solveMaster", "failed to solve the master problem", 0);
		return 1;
	}

	getCallbackPrimal(cbdata, wherefrom, args->cell->candidX, args->prob[0]->num->cols);

#if defined(UserMIPcutsActive)
	/* Use GMI and MIR cutting planes to solve the SD-optimized problem */
	args->cell->MIPFlag = false;
	if (solveIntCell(args->stoc, args->prob, args->cell)) {
		errMsg("algorithm", "algo", "failed to solve the cell using GMI and MIR algorithm", 0);
		return 1;
	}
#endif // defined(UserMIPcutsActive)


	/*Update the number of cuts*/
	args->MIPcuts += getNumRows(args->cell->master->lp) - args->OPTcuts;

	/* Invoke the main loop of the L-shaped method */
	while (args->cell->kii < config.MAX_ITER_CLBK) {
		tic = clock();
		args->OPTcuts++;

		printf("%i - %i", args->cell->kii, config.MAX_ITER_CLBK);

		if (mainloopSDCell_callback(args->stoc, args->prob, args->cell, &breakLoop, observ) ) {
			errMsg("Callback", "usersolve", "failed to solve Benders cell for the node problem", 0);
			return 1;
		}


		args->cell->time.masterAccumTime += args->cell->time.masterIter; args->cell->time.subprobAccumTime += args->cell->time.subprobIter;
		args->cell->time.argmaxAccumTime += args->cell->time.argmaxIter; args->cell->time.optTestAccumTime += args->cell->time.optTestIter;
		args->cell->time.masterIter = args->cell->time.subprobIter = args->cell->time.optTestIter = 0.0;
		args->cell->time.argmaxIter = args->cell->time.optTestIter = 0.0;
		args->cell->time.iterTime = ((double) clock() - tic)/CLOCKS_PER_SEC; args->cell->time.iterAccumTime += args->cell->time.iterTime;
		printf("#"); fflush(stdout);

		if ( breakLoop )
			break;
	}

	if (args->call < config.NodeNum)
	{
		args->cell->nodeSol[args->call].sol = duplicVector(args->cell->candidX, args->prob[0]->num->cols);
		args->cell->nodeSol[args->call].sol_size = args->prob[0]->num->cols;
		args->cell->nodeSol[args->call].nodeNum = args->call;
		args->cell->nodeSol[args->call].LB = args->cell->candidEst;
		args->cell->nodeSol[args->call].piM = duplicVector(args->cell->piM, args->cell->master->mar);
		args->cell->nodeSol[args->call].djM = duplicVector(args->cell->djM, args->cell->master->mar);
		args->cell->nodeSol[args->call].mar = args->cell->master->mar;
		args->cell->nodeSol[args->call].isInt = isInteger(args->cell->candidX, args->prob[0]->num->cols,
			0, args->prob[0]->num->cols - 1, config.INT_TOLERANCE);
	}

	args->cell->master->lp = NULL;
	args->cell->master->lp = nodelp;
	mem_free(candidCuts);

	return 0;
} /* END usersolve */


