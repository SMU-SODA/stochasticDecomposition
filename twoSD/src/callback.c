/*
 * callback.c
 *
 *  Created on: Feb 24, 2021
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"

int phase_one_analysis(stocType *stoc, probType **prob, cellType *cell) {

	int status = 0;
	printf("\n--------------------------------Phase 1 summary------------------------------------------\n");
	//0a - (LB on xLP is cell->incumbEst)
	//0b - xLP is phase1 solution -> evaluate xLP (UB)
	/* Evaluate the compromise solution */
	printf("LB estimate: %0.4f\n", cell->incumbEst);
	evaluate(NULL,NULL, stoc, prob, cell->subprob, cell->candidX);


	//1a - Change the master to MILP no callback solve to optimality - > xIP (LB)
	//1b - check the optimality of (xIP) using pretest
	//1c - Evaluate the solution (xIP) using evaluate(sFile, stoc, prob, cell->subprob, cell->incumbX); -> (UB)
	/* Turn the clone problem to LP */
	QPtoLP(stoc, prob, cell, 0);

	/* Launch the solver to solve the MIP master problem */
	if (solveProblem(cell->master->lp, cell->master->name, cell->master->type, &status)) {
		errMsg("algorithm", "algo-after-phase1", "failed to solve the master problem", 0);
		return 1;
	}
	/* Get the most recent optimal solution to master program */
	if (getPrimal(cell->master->lp, cell->candidX, prob[0]->num->cols)) {
		errMsg("algorithm", "solveMaster", "failed to obtain the primal solution for master", 0);
		return 1;
	}

	/* Find the highest cut at the candidate solution. where cut_height = alpha - beta . x */ //
	cell->candidEst = vXvSparse(cell->candidX, prob[0]->dBar) + maxCutHeight(cell->activeCuts, cell->sampleSize, cell->candidX, prob[0]->num->cols, prob[0]->lb);
	printf("\n\nLP estimate: %0.4f\n", cell->candidEst);
	printVector(cell->candidX, prob[0]->num->cols, NULL);
	evaluate(NULL,NULL, stoc, prob, cell->subprob, cell->candidX);

	/* Turn the clone problem to MILP */
	LPtoMILP(stoc, prob, cell);
	/* Launch the solver to solve the MIP master problem */
	if (solveProblem(cell->master->lp, cell->master->name, cell->master->type, &status)) {
		errMsg("algorithm", "algo-after-phase1", "failed to solve the master problem", 0);
		return 1;
	}
	/* Get the most recent optimal solution to master program */
	if (getPrimal(cell->master->lp, cell->candidX, prob[0]->num->cols)) {
		errMsg("algorithm", "solveMaster", "failed to obtain the primal solution for master", 0);
		return 1;
	}

	/* Find the highest cut at the candidate solution. where cut_height = alpha - beta . x */ //
	cell->candidEst = vXvSparse(cell->candidX, prob[0]->dBar) + maxCutHeight(cell->activeCuts, cell->sampleSize, cell->candidX, prob[0]->num->cols, prob[0]->lb);
	printf("\n\nMILP estimate: %0.4f\n", cell->candidEst);
	printVector(cell->candidX, prob[0]->num->cols, NULL);
	evaluate(NULL,NULL, stoc, prob, cell->subprob, cell->candidX);

	//summary for the second one

	printf("----------------------------------------------------------------------------------------------------\n");

	return 0;
}


/*
  Turn the QP master problem to LP
  Siavash Tabrizian May 20
*/
int QPtoLP(stocType *stoc, probType **prob, cellType *cell, int toMIP) {

	cString lu; cString uu;
	iVector indices;
	int numCols = prob[0]->num->cols;
	int status = 0;

	if (!(indices = (iVector)arr_alloc(numCols, int)))
		errMsg("allocation", "bendersCallback", "indices", 0);
	for (int c = 0; c < numCols; c++)
		indices[c] = c;

	if (toMIP == 0)
	{

		/*********00. Set sigma to 0 *********/
		if (changeQPproximal(cell->master->lp, prob[0]->num->cols, 0.0)) {
			errMsg("algorithm", "SolveIntCell", "failed to change the proximal term", 0);
			return 1;
		}

		/*********01. Update the right hand sides of master *********/
		if (revchangeQPrhs(prob[0], cell, cell->incumbX)) {
			errMsg("algorithm", "SolveIntCell", "failed to change the right-hand side to convert the problem into QP", 0);
			return 1;
		}

		/*********02. Change LP solver to primal simplex *********/
		if (changeProbType(cell->master->lp, PROB_LP))
		{
			errMsg("algorithm", "SolveIntCell", "failed to set primal algorithm for LP", 0);
			return 1;
		}
		if (changeLPSolverType(PROB_LP))
		{
			errMsg("algorithm", "SolveIntCell", "failed to set primal algorithm for LP", 0);
			return 1;
		}

	}
	else
	{
		/*********00. Set sigma to 0 *********/
		//if (changeQPproximal(cell->master->lp, prob[0]->num->cols, 0.0)) {
		//	errMsg("algorithm", "SolveIntCell", "failed to change the proximal term", 0);
		//	return 1;
		//}

		/*********01. Update the right hand sides of master *********/
		if (revchangeQPrhs(prob[0], cell, cell->incumbX)) {
			errMsg("algorithm", "SolveIntCell", "failed to change the right-hand side to convert the problem into QP", 0);
			return 1;
		}

		cell->master->type = PROB_MILP;

		if (changeCtype(cell->master->lp, numCols, indices, cell->master->ctype)) {
			errMsg("solver", "bendersCallback", "failed to change column type", 0);
			return 1;
		}

		if (changeProbType(cell->master->lp, PROB_MILP)) {
			errMsg("Problem Setup", "bendersCallback", "master", 0);
			return 1;
		}

		/*********02. Change LP solver to B&B *********/
		if (changeLPSolverType(PROB_MILP))
		{
			errMsg("algorithm", "SolveIntCell", "failed to set primal algorithm for LP", 0);
			return 1;
		}

	}

	///* Change the column types to integer/binary to change the problem type to SMIP */
	if (!(lu = (cString)arr_alloc(numCols, int)))
		errMsg("allocation", "algo-copymaster", "lu", 0);
	if (!(uu = (cString)arr_alloc(numCols, int)))
		errMsg("allocation", "algo-copymaster", "uu", 0);


	for (int c = 0; c < numCols; c++)
	{
		lu[c] = 'L';
		uu[c] = 'U';
	}
	status = changeBDS(cell->master->lp, numCols, indices, lu, prob[0]->sp->bdl);
	if (status) {
		solverErrmsg(status);
		errMsg("solver", "algo-copymaster", "failed to change column lowerbound", 0);
		return 1;
	}
	status = changeBDS(cell->master->lp, numCols, indices, uu, prob[0]->sp->bdu);
	if (status) {
		solverErrmsg(status);
		errMsg("solver", "algo-copymaster", "failed to change column upperbound", 0);
		return 1;
	}


#if defined(LPMIP_PRINT)
	if (toMIP == 0)
	{
		writeProblem(cell->master->lp, "finalMaster_QP2LP.lp");
	}
	else
	{
		writeProblem(cell->master->lp, "finalMaster_QP2MILP.lp");
	}
#endif

	mem_free(indices);

	return 0;
}

/*
Turn the LP master problem to MILP
Siavash Tabrizian July 20
*/
int LPtoMILP(stocType *stoc, probType **prob, cellType *cell) {
	iVector indices;
	int numCols = prob[0]->num->cols;

	if (!(indices = (iVector)arr_alloc(numCols, int)))
		errMsg("allocation", "bendersCallback", "indices", 0);
	for (int c = 0; c < numCols; c++)
		indices[c] = c;

	cell->master->type = PROB_MILP;

	if (changeCtype(cell->master->lp, numCols, indices, cell->master->ctype)) {
		errMsg("solver", "bendersCallback", "failed to change column type", 0);
		return 1;
	}

	if (changeProbType(cell->master->lp, PROB_MILP)) {
		errMsg("Problem Setup", "bendersCallback", "master", 0);
		return 1;
	}

	/*********02. Change LP solver to B&B *********/
	if (changeLPSolverType(PROB_MILP))
	{
		errMsg("algorithm", "SolveIntCell", "failed to set primal algorithm for LP", 0);
		return 1;
	}


#if defined(LPMIP_PRINT)
	writeProblem(cell->master->lp, "finalMaster_LP2MILP.lp");
#endif

	mem_free(indices);
	return 0;
}//END LPtoMILP();

void printNodeInfo(nodeInfo    *nodeSol, int Nodecnt)
{
	printf("\n --------------- Printing the node solutions Info: \n");
	printf("\n number of nodes: %i \n", Nodecnt);
	for (int i = 0; i < Nodecnt; i++)
	{
		printf("\n ----- ------- ----- \n");
		printf("\n call: %i \n", nodeSol[i].nodeNum);
		printVector(nodeSol[i].sol, nodeSol[i].sol_size, NULL);
		printf("\n Estimate: %0.4f", nodeSol[i].LB);
		if (nodeSol[i].isInt)
		{
			printf("\n integer solution!\n");
		}
		printf("\n dual of the master for this node: \n");
		printVector(nodeSol[i].piM, nodeSol[i].mar, NULL);
		printf("\n ----- ------- ----- \n");
	}
	printf("\n ---------------------------------------------------\n");
}//END printNodeInfo()
