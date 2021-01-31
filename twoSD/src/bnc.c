/*
* bnc.c
*
*  Created on: Nov 30, 2020
*      Author: Siavash Tabrizian 
* Institution: Southern Methodist University
*
* Please send you comments or bug report to stabrizian (at) smu (dot) edu
*
* This is part of the SD-integer project which handels the branch and bound tree
*
*/

#include "twoSD.h"


extern configType config;
#define maxdnodes   1000
#define useDNODE 
#define useINODE



#undef writeprob

int currKey;
int currDepth;
double GlobeUB;               // Global upper bound
oneProblem      *original;    // Info of the original problem 
struct BnCnodeType *bestNode; // best node that is found so far
struct BnCnodeType **nodearr; // array of deactivated leaf nodes
struct BnCnodeType **inodearr;// array of integer feasible leaf nodes
int dnodes;                   // number of deactivated nodes 
int inodes;                   // number of integer feasible nodes 
int maxcut;
double meanVal;               // Global lower bound 
#define printBest
#define printBranch
#undef depthtest
#undef writemaster


int sumintVec(iVector a, int len)
{
	int n;
	int output = 0.0;

	for (n = 0; n < len; n++)
		output += a[n];

	return output;
}

/* check if an integer is in an iVector */
bool isInVec(iVector vec, int len, int val)
{
	int i;

	for (i = 0; i < len; i++)
		if (vec[i] == val)
		{
			return true;
		}

	return false;
}

struct BnCnodeType *newrootNode(int numVar, double LB, double UB, oneProblem * orig)
{
	int i; int j;

	struct BnCnodeType *temp = (struct BnCnodeType *)malloc(sizeof(struct BnCnodeType));
	temp->key = 0;
	temp->parentkey = -1;
	temp->depth = 0;
	temp->isleft = false;
	temp->fracVal = -1;
	temp->varId = -1;
	temp->ishrstic = false;
	temp->left = temp->right = NULL;
	temp->nextnode = temp->prevnode = NULL;
	temp->numVar = numVar;
	temp->LB = LB;
	temp->UB = UB;
	temp->parentnumSamp = 0;
	temp->numSamp = 0;
	temp->partightPi = 0;
	temp->parLambdasize = 0;
	temp->tightPi = 0;
	temp->Lambdasize = 0;
	temp->parparinit = 0;
	temp->fracPi = 0.0;
	temp->isActive = true;
	temp->isfathomed = false;
	temp->parobjVal = INFINITY;
	temp->isSPopt = true;
	if (!(temp->disjncs = (iVector)arr_alloc(temp->numVar, double)))
		errMsg("allocation", "newNode", "temp->disjncs", 0);
	if (!(temp->disjncsVal = (dVector *)arr_alloc(temp->numVar, dVector)))
		errMsg("allocation", "newNode", "temp->disjncs", 0);
	for (i = 0; i < numVar; i++)
	{
		if (!(temp->disjncsVal[i] = (dVector)arr_alloc(2, double)))
			errMsg("allocation", "newNode", "temp->disjncs", 0);
	}
	if (!(temp->vars = (dVector)arr_alloc(temp->numVar + 1, double)))
		errMsg("allocation", "newNode", "temp->vars", 0);
	if (config.Pi_EVAL_FLAG == 1)
	{
		if (!(temp->IncumbiStar = (iVector)arr_alloc(maxcut*config.MAX_ITER, int)))
			errMsg("allocation", "newNode", "temp->vars", 0);
		if (!(temp->ParIncumbiStar = (iVector)arr_alloc(maxcut*config.MAX_ITER, int)))
			errMsg("allocation", "newNode", "temp->vars", 0);
		for (i = 0; i < maxcut*config.MAX_ITER; i++)
		{
			temp->IncumbiStar[i] = -1;
			temp->ParIncumbiStar[i] = -1;
		}
	}
	for (int v = 0; v < temp->numVar + 1; v++)
		temp->vars[v] = 0.0;
	for (i = 0; i < temp->numVar; i++)
	{
		temp->disjncsVal[i][0] = orig->bdl[i];
		temp->disjncsVal[i][1] = orig->bdu[i];
		temp->disjncs[i] = 0;
	}

	return temp;

}//End newrootNode()


struct BnCnodeType *newNode(int key, struct BnCnodeType * parent, double fracVal, int varId, bool isleft)
{
	int i; int j;
	if (fabs(fracVal - round(fracVal)) == 0)
	{
		if (round(fracVal) == parent->disjncsVal[varId][0]) {
			fracVal += 0.001;
		}
		else if (round(fracVal) == parent->disjncsVal[varId][1]) {
			fracVal -= 0.001;
		}
		else {
			fracVal += 0.001;
		}
	}
	struct BnCnodeType *temp = (struct BnCnodeType *)malloc(sizeof(struct BnCnodeType));
	temp->key = key;
	temp->parentkey = parent->key;
	temp->depth = parent->depth + 1;
	temp->isleft = isleft;
	temp->fracVal = fracVal;
	temp->varId = varId;
	temp->ishrstic = false;
	temp->left = temp->right = NULL;
	temp->nextnode = temp->prevnode = NULL;
	temp->numVar = parent->numVar;
	temp->stInt = parent->stInt;
	temp->edInt = parent->edInt;
	temp->LB = parent->LB;
	temp->UB = parent->UB;
	temp->parentnumSamp = parent->numSamp;
	temp->numSamp = 0;
	temp->partightPi = parent->tightPi;
	temp->parLambdasize = parent->Lambdasize;
	temp->tightPi = 0;
	temp->Lambdasize = 0;
	temp->parparinit = parent->parLambdasize;
	temp->fracPi = 0.0;
	temp->isActive = true;
	temp->isfathomed = false;
	temp->parobjVal = parent->objval;
	temp->isSPopt = true;
	if (!(temp->disjncs = (iVector)arr_alloc(temp->numVar, int)))
		errMsg("allocation", "newNode", "temp->disjncs", 0);
	if (!(temp->disjncsVal = (dVector *)arr_alloc(temp->numVar, dVector)))
		errMsg("allocation", "newNode", "temp->disjncs", 0);
	if (config.Pi_EVAL_FLAG == 1)
	{
		if (!(temp->IncumbiStar = (iVector)arr_alloc(maxcut*config.MAX_ITER, int)))
			errMsg("allocation", "newNode", "temp->vars", 0);
		if (!(temp->ParIncumbiStar = (iVector)arr_alloc(maxcut*config.MAX_ITER, int)))
			errMsg("allocation", "newNode", "temp->vars", 0);
		for (i = 0; i < maxcut*config.MAX_ITER; i++)
		{
			temp->IncumbiStar[i] = -1;
			temp->ParIncumbiStar[i] = parent->IncumbiStar[i];
		}
	}

#if defined(BNC_CHECK)
	printf("\nfrac val: %0.4f\n",fracVal);
#endif // defined(BNC_CHECK)

	for (i = 0; i < parent->numVar; i++)
	{
		if (!(temp->disjncsVal[i] = (dVector)arr_alloc(2, double)))
			errMsg("allocation", "newNode", "temp->disjncs", 0);
		if (i == varId)
		{
			if (temp->isleft == true)
			{
				if (config.BRN_STR == 0)
				{
					temp->disjncs[i] = 1;
					temp->disjncsVal[i][0] = parent->disjncsVal[i][0];
					temp->disjncsVal[i][1] = floor(fracVal);
#if defined(BNC_CHECK)
					printf("\nis left: %d - var id: %d - ub: %0.4f\n", temp->isleft, i, floor(fracVal));
#endif // defined(BNC_CHECK)
				}
				else if (config.BRN_STR == 1)
				{
					if (abs(fracVal - ceil(fracVal)) > 0.5)
					{
						temp->disjncs[i] = 1;
						temp->disjncsVal[i][0] = ceil(fracVal);
						temp->disjncsVal[i][1] = parent->disjncsVal[i][1];
					}
					else
					{
						temp->disjncs[i] = 1;
						temp->disjncsVal[i][0] = parent->disjncsVal[i][0];
						temp->disjncsVal[i][1] = floor(fracVal);
					}
				}
			}
			else
			{
				if (config.BRN_STR == 0)
				{
					temp->disjncs[i] = 1;
					temp->disjncsVal[i][0] = ceil(fracVal);
					temp->disjncsVal[i][1] = parent->disjncsVal[i][1];
#if defined(BNC_CHECK)
					printf("\nis left: %d - var id: %d - lb: %0.4f\n", temp->isleft, i, ceil(fracVal));
#endif // defined(BNC_CHECK)
				}
				else if (config.BRN_STR == 1)
				{
					if (abs(fracVal - ceil(fracVal)) <= 0.5)
					{
						temp->disjncs[i] = 1;
						temp->disjncsVal[i][0] = ceil(fracVal);
						temp->disjncsVal[i][1] = parent->disjncsVal[i][1];
					}
					else
					{
						temp->disjncs[i] = 1;
						temp->disjncsVal[i][0] = parent->disjncsVal[i][0];
						temp->disjncsVal[i][1] = floor(fracVal);
					}
				}
			}
		}
		else
		{
			temp->disjncs[i] = parent->disjncs[i];
			temp->disjncsVal[i][0] = parent->disjncsVal[i][0];
			temp->disjncsVal[i][1] = parent->disjncsVal[i][1];
		}
	}
	if (!(temp->vars = (dVector)arr_alloc(temp->numVar + 1, double)))
		errMsg("allocation", "newNode", "temp->vars", 0);

	for (int v = 0; v < temp->numVar + 1; v++)
		temp->vars[v] = parent->vars[v];

	return temp;

}//End newNode()

struct BnCnodeType *copyNode(struct BnCnodeType *node, double thresh)
{
	int i; int j;
	struct BnCnodeType *temp = (struct BnCnodeType *)malloc(sizeof(struct BnCnodeType));

	temp->key = node->key;
	temp->parentkey = node->parentkey;
	temp->depth = node->depth;
	temp->isleft = node->isleft;
	temp->fracVal = node->fracVal;
	temp->ishrstic = true;
	temp->varId = node->varId;
	temp->left = temp->right = NULL;
	temp->nextnode = temp->prevnode = NULL;
	temp->numVar = node->numVar;
	temp->stInt = node->stInt;
	temp->edInt = node->edInt;
	temp->LB = node->LB;
	temp->UB = node->UB;
	temp->parentnumSamp = node->parentnumSamp;
	temp->numSamp = node->numSamp;
	temp->partightPi = node->partightPi;
	temp->parLambdasize = node->parLambdasize;
	temp->tightPi = node->tightPi;
	temp->Lambdasize = node->Lambdasize;
	temp->parparinit = node->parparinit;
	temp->fracPi = node->fracPi;
	temp->isActive = true;
	temp->isfathomed = false;
	temp->parobjVal = node->parobjVal;
	temp->isSPopt = node->isSPopt;
	if (!(temp->disjncs = (iVector)arr_alloc(temp->numVar, int)))
		errMsg("allocation", "newNode", "temp->disjncs", 0);
	if (!(temp->disjncsVal = (dVector *)arr_alloc(temp->numVar, dVector)))
		errMsg("allocation", "newNode", "temp->disjncs", 0);
	if (!(temp->vars = (dVector)arr_alloc(temp->numVar + 1, double)))
		errMsg("allocation", "newNode", "temp->vars", 0);
	if (config.Pi_EVAL_FLAG == 1)
	{
		if (!(temp->IncumbiStar = (iVector)arr_alloc(maxcut*config.MAX_ITER, int)))
			errMsg("allocation", "newNode", "temp->vars", 0);
		if (!(temp->ParIncumbiStar = (iVector)arr_alloc(maxcut*config.MAX_ITER, int)))
			errMsg("allocation", "newNode", "temp->vars", 0);
		for (i = 0; i < maxcut*config.MAX_ITER; i++)
		{
			temp->ParIncumbiStar[i] = node->ParIncumbiStar[i];
			temp->IncumbiStar[i] = node->IncumbiStar[i];
		}
	}
	for (i = 0; i < node->numVar; i++)
	{
		if (!(temp->disjncsVal[i] = (dVector)arr_alloc(2, double)))
			errMsg("allocation", "newNode", "temp->disjncs", 0);
		double val = 0.0;
		if (node->vars[i] > thresh) val = ceil(node->vars[i]); else val = floor(node->vars[i]);
		temp->disjncs[i] = 1;
		temp->disjncsVal[i][0] = val;
		temp->disjncsVal[i][1] = val;
		temp->vars[i] = val;
	}
	
	return temp;
}//End copyNode()


 /* after creating the node problem (lp) we can impose the disjunctions as new bounds on
 variables */
int addBnCDisjnct(cellType *cell, dVector  *disjncsVal, int numCols, struct BnCnodeType * node)
{
	int 	status = 0, cnt;
	dVector	lbounds, ubounds;

	if (!(lbounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "addDisjnct", "lbounds", 0);

	if (!(ubounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "addDisjnct", "ubounds", 0);

	/* Change the Upper Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		ubounds[cnt] = disjncsVal[cnt][1];
		cell->master->bdu[cnt] = disjncsVal[cnt][1];
	}

	/* Change the Lower Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		lbounds[cnt] = disjncsVal[cnt][0];
		cell->master->bdl[cnt] = disjncsVal[cnt][0];
	}

	if (!node->ishrstic && node->depth != 0)
	{
		if (node->isleft)
		{
			node->vars[node->varId + 1] = disjncsVal[node->varId][1];
		}
		else {
			node->vars[node->varId + 1] = disjncsVal[node->varId][0];
		}
	}

#if defined(BNC_CHECK)
	printf("\nlower bounds:\n");
	printVector(cell->master->bdl, numCols, NULL);
	printf("\nupper bounds:\n");
	printVector(cell->master->bdu, numCols, NULL);
#endif // defined(BNC_CHECK)

	
	if (changeQPbds(cell->master->lp, numCols, lbounds, ubounds, node->vars, 0)) {
		errMsg("algorithm", "algoIntSD", "failed to change the bounds to convert the problem into QP", 0);
		return 1;
	}

#if defined(BNC_CHECK)
	printf("\nafter changeQPbds");
	if (getLb(cell->master->lp, 0, numCols, lbounds)) {
		errMsg("bnc", "addBnCDisjnct", "failed to get lb", 0);
		return 1;
    }
	if (getUb(cell->master->lp, 0, numCols, ubounds)) {
		errMsg("bnc", "addBnCDisjnct", "failed to get lb", 0);
		return 1;
	}
	printf("\nlower bounds:\n");
	printVector(lbounds, numCols, NULL);
	printf("\nupper bounds:\n");
	printVector(ubounds, numCols, NULL);
#endif // defined(BNC_CHECK)

	mem_free(lbounds); 
	mem_free(ubounds);

	return 0;
}

// Truncate the var from SD based on the lower bounds and upper bounds of the original problem
void truncate(dVector var, dVector lb, dVector ub, int cnt)
{
	for (int v = 0; v < cnt; v++)
	{
		if (var[v+1] < lb[v]) var[v+1] = lb[v];
		if (var[v+1] > ub[v]) var[v+1] = ub[v];
	}
}

 // Subroutine for Solving the node problem given the lp pointer and 
 // node information
double solveNode(stocType *stoc, probType **prob, cellType *cell, struct BnCnodeType *node, cString pname)
{
	int 	status;
	cell->ki = 0;
	
#if defined(BNC_CHECK)
	printf("\nbefore SD: %-10s%-5d%-10s%-5d%-5s%-5d%-10s%-5d%-10s%-5d","node id:",
		node->key,"node depth:",node->depth,"numVar:",node->numVar,"node->edInt:",node->edInt,"master->mac:",cell->master->mac);
	printf("\ndisjuction:\n");
	printIntvec(node->disjncs, node->numVar-1, NULL);
#endif // defined(BNC_CHECK)

#if defined(clean_master)
	if (node->depth > 0)
	{
		/* Get the total number of rows */
		int row_num = getNumRows(cell->master->lp);

#if defined(writemaster)
		writeProblem(cell->master->lp, "master_test_beforeclean.lp");
#endif // defined(writemaster)


		if (row_num > prob[0]->num->rows)
		{
			/* Get rid of the indexed cut on the solver */
			/* oneProblem structures and solver elements */
			if (prob[0]->num->rows < row_num)
			{
				for (int cnt = row_num - 1; cnt >= prob[0]->num->rows; cnt--)
					if (removeRow(cell->master->lp, cnt, cnt)) {
						printf("row Num %d - tot rows %d - orig rows %d", cnt, row_num, prob[0]->num->rows);
						errMsg("solver", "cleanCellType", "failed to remove a row from master problem", 0);
						return 1;
					}
			}
		}

		/* deactivate the current cuts */
		for (int c = 0; c < cell->cuts->cnt; c++)
		{
			cell->cuts->vals[c]->isAvctive = false;
		}

		cell->master->mar = prob[0]->num->rows-1;

#if defined(writemaster)
		writeProblem(cell->master->lp, "master_test_afterclean.lp");
#endif // defined(writemaster)

	}
#endif // defined(clean_master)


	if (true)
	{

		if (addBnCDisjnct(cell, node->disjncsVal, node->edInt, node))
			errMsg("addDisjnct", "solveNode", "adding disjunctions are failed", 0);

		////Initializing candidX, candidEst, incumbEst and IncumbX
		truncate(node->vars, prob[0]->sp->bdl, prob[0]->sp->bdu, node->numVar);
		copyVector(node->vars, cell->incumbX, node->numVar, true);
		copyVector(cell->incumbX, cell->candidX, node->numVar, true);
		cell->candidEst = vXvSparse(cell->candidX, prob[0]->dBar) + maxCutHeight(cell->cuts, cell->sampleSize, cell->candidX, prob[0]->num->cols, prob[0]->lb);
		cell->incumbEst = cell->candidEst;
		node->objval = cell->incumbEst;

		
		double testval = maxCutHeight(cell->cuts, cell->sampleSize, cell->candidX, prob[0]->num->cols, prob[0]->lb);
		// Setup the QP master problem 
		if (changeQPproximal(cell->master->lp, node->edInt, cell->quadScalar)) {
			errMsg("algorithm", "algoIntSD", "failed to change the proximal term", 0);
			return 1;
		}

		if (changeQPrhs(prob[0], cell, node->vars)) {
			errMsg("algorithm", "algoIntSD", "failed to change the right-hand side to convert the problem into QP", 0);
			return 1;
		}


	}

#if defined(BNC_CHECK)
	printf("\ninit var:\n");
	printVector(node->vars, node->numVar, NULL);
	printf("\ninit incumbX:\n");
	printVector(cell->incumbX, node->numVar, NULL);
#endif // defined(BNC_CHECK)
	
	int k_before = cell->k;
	if (node->ishrstic || (cell->k < config.MAX_ITER && (node->prevnode == NULL || (node->objval < GlobeUB && !isInteger(node->vars, node->edInt, 0, node->edInt + 1, config.TOLERANCE)))))
	{
#if defined(writemaster)
		writeProblem(cell->master->lp, "master_test_beforesolve.lp");
#endif // defined(writemaster)

		/* Use two-stage stochastic decomposition algorithm to solve the problem */
		if (solveCell(stoc, prob, cell)) {
			return 1;
		}
		int k_after = cell->k;
		int k_diff = k_after - k_before;

#if defined(BNC_CHECK)
		printf("\nafter SD:%-10s%-10.4f%-10s%-10.4f", "incumbEst:", cell->incumbEst, "candidEst:", cell->candidEst);
		printf("\nincumbX SD:\n");
		printVector(cell->incumbX, cell->master->mac, NULL);
		printf("\ncandidX SD:\n");
		printVector(cell->candidX, cell->master->mac, NULL);
#endif // defined(BNC_CHECK)

		copyVector(cell->incumbX, node->vars, node->numVar, true);
		truncate(node->vars, prob[0]->sp->bdl, prob[0]->sp->bdu, node->numVar);

#if defined(writemaster)
		writeProblem(cell->master->lp, "master_test_aftersolve.lp");
#endif // defined(writemaster)

		if (cell->incumbEst <= node->parobjVal || cell->incumbEst <= meanVal)
		{
			node->objval = meanVal - fabs(cell->incumbEst <= meanVal);
			node->isSPopt = false;
		}
		else
		{
			node->objval = cell->incumbEst;
		}
		cell->optFlag = false;

#if defined(BNC_CHECK)
		cString nullString = NULL;
		printf("\nafter SD var:\n");
		printVector(node->vars, node->numVar, NULL);
		/*nullString = (cString)arr_alloc(2 * BLOCKSIZE, char);
		printf(">>>press a key to continue...");
		scanf("%s", (*nullString));*/
#endif // defined(BNC_CHECK)

	}
	else
	{
		truncate(node->vars, prob[0]->sp->bdl, prob[0]->sp->bdu, node->numVar);
		if (cell->incumbEst <= node->parobjVal || cell->incumbEst <= meanVal)
		{
			node->objval = meanVal - fabs(cell->incumbEst <= meanVal);
			node->isSPopt = false;
		}
		else
		{
			node->objval = cell->incumbEst;
		}
		printf("SD output: iters:%-4d - dnodes:%-3d - sigma size:%-7d - lambda size:%-7d - omega size:%-7d\n", cell->ki, dnodes + 1,cell->sigma->cnt,cell->lambda->cnt, cell->omega->cnt);

		return 0;
	}

	if (cell->ki > 600) printf("\n");;
	printf("SD output: iters:%-4d - dnodes:%-3d - sigma size:%-7d - lambda size:%-7d - omega size:%-7d\n", cell->ki, dnodes + 1, cell->sigma->cnt, cell->lambda->cnt, cell->omega->cnt);

	return 0;
}


// Insert a node based on its key to the tree 
struct BnCnodeType *insertNode(struct BnCnodeType *node, struct BnCnodeType *activenode) {
	// Return an error message if the inserted node is NULL
	if (node == NULL)
	{
		errMsg("NULL node", "insert", "node is NULL", 0);
	}

	// Traverse to place the node - insert it to the queue
	node->prevnode = activenode;
	activenode->nextnode = node;

	return node;
}

// return the inorder successor of the input node: this is just based on the order 
// of the node in the tree not lower bound or upper bound
struct BnCnodeType *successorNode(struct BnCnodeType *node) {
	struct BnCnodeType *current = node;

	// Find the leftmost leaf
	while (current && current->left != NULL)
		current = current->left;

	return current;
}

int freeNodes(struct BnCnodeType *root)
{
	// Return 1 when the tree is empty
	if (root == NULL) return 1;

	struct BnCnodeType *next = root;

	while (next)
	{
		struct BnCnodeType *node = next;
		next = node->nextnode;
		if (node) {
			if (node->disjncs) mem_free(node->disjncs);
			if (node->disjncsVal)  mem_free(node->disjncsVal);
			if (node->vars)  mem_free(node->vars);
			//if (node->duals)  mem_free(node->duals);
			//if (node->prevnode) node->prevnode->nextnode = NULL;
			//if (node->nextnode) node->nextnode->prevnode = NULL;
			mem_free(node);
		}
	}


	return 0;

}//End freeNode()

int freeNode(struct BnCnodeType *node)
{
	// Return 1 when the tree is empty
	if (node == NULL) return 1;

	if (node->prevnode) node->prevnode->nextnode = NULL;
	if (node->nextnode) node->nextnode->prevnode = NULL;
	node->prevnode = NULL; node->nextnode = NULL;
	if (node->disjncs) mem_free(node->disjncs);
	if (node->disjncsVal)  mem_free(node->disjncsVal);
	if (node->vars)  mem_free(node->vars);
	//if (node->duals)  mem_free(node->duals);
	mem_free(node);


	return 0;

}//End freeNode()


 /* given a node problem which is solved decide which variable should be selected
 for branching  */
int branchVar(struct BnCnodeType *node, int strategy)
{

	if (strategy == 0) // strategy based on the least index of the fractional variable 
	{
		return node->depth;
	}
	else if (strategy == 1)// strategy based on the larger fractional value
	{
		double larger = -1;  /* maximum fractional value */
		int maxidx = node->depth;

		for (int i = 0; i < node->edInt; i++)
		{
			if (node->disjncsVal[i][0]!= node->disjncsVal[i][1] && larger < node->vars[i])
			{
				larger = node->vars[i];
				maxidx = i;
			}
		}

		return maxidx;
	}
	else if (strategy == 2)// strategy based on the smaller fractional value
	{
		double smaller = INFINITY;  /* minimum fractional value */
		int minidx = node->depth;

		for (int i = 0; i < node->edInt; i++)
		{
			if (node->disjncsVal[i][0] != node->disjncsVal[i][1] && smaller > node->vars[i])
			{
				smaller = node->vars[i];
				minidx = i;
			}
		}

		return minidx;
	}
	else if (strategy == 3)// strategy based on the value closer to 0.5
	{
		double avg = INFINITY;  /* avg fractional value */
		int idx = node->depth;

		for (int i = 0; i < node->edInt; i++)
		{
			if (node->disjncsVal[i][0] != node->disjncsVal[i][1] && abs(0.5 - node->vars[i]) < avg)
			{
				avg = abs(0.5 - node->vars[i]);
				idx = i;
			}
		}

		return idx;
	}

}

/* branch and bound algorithm */
int branchbound(stocType *stoc, probType **prob, cellType *cell, double LB, double UB)
{

	int i, j;
	int nodecnt = 0;
	int maxcut = config.CUT_MULT * cell->master->mac + 3;
	cell->basis->incumPicnt = maxcut;
	cell->basis->basisEval = config.Pi_EVAL_FLAG;

#if defined(useDNODE)
	if (!(nodearr = (struct BnCnodeType **)arr_alloc(maxdnodes, struct BnCnodeType *)))
		errMsg("allocation", "branchbound", "nodearr", 0);
#endif // defined(useDNODE)
#if defined(useINODE)
	if (!(inodearr = (struct BnCnodeType **)arr_alloc(maxdnodes, struct BnCnodeType *)))
		errMsg("allocation", "branchbound", "inodearr", 0);
#endif // defined(useDNODE)
	dnodes = -1;
	inodes = -1;

	/* set of LB and UB */
	GlobeUB = UB;

	original = prob[0]->sp;
	struct BnCnodeType * rootNode = NULL;  // root node
	struct BnCnodeType * activeNode = NULL;// active node in the tree 
	struct BnCnodeType * currentNode = NULL;// current node that is investigated

#if defined(printBranch)
	printLine();
	printf("Starting the BnC branch and bound procedure...\n");
	printf("Github page: https://github.com/siavashtab \n");
	printf("Email: stabrizian@smu.edu\n");
	printLine();
	printf("\n\n");
	printLongLine();
	printf("%-10s%-10s%-10s%-10s%-10s%-12s%-12s%-12s%-12s%-12s\n", "node id", "parent", "depth", "k", "\|x\|","fval", "UB", "feasible", "integer","Lamfrac");
	printLongLine();
#endif // defined(printBranch)

	int totdepth = original->mac;
	double totnodes = pow(2, totdepth);

	rootNode = newrootNode(original->mac, LB, UB, original);
	rootNode->parobjVal = cell->meanVal;
	meanVal = cell->meanVal;

	for (i = 0; i < rootNode->numVar; i++)
	{
		if (original->ctype[i] == 'B' || original->ctype[i] == 'I') {
			rootNode->stInt = i;
			break;
		}
	}
	rootNode->edInt = 1;
	for (i = rootNode->stInt + 1; i < rootNode->numVar; i++)
	{
		if (original->ctype[i] != 'B' && original->ctype[i] != 'I') {
			rootNode->edInt = i;
			break;
		}
		else
		{
			rootNode->edInt++;
		}
	}
	for (i = 0; i < prob[0]->num->cols; i++)
		rootNode->vars[i] = cell->incumbX[i+1];
	activeNode = rootNode;
	currentNode = rootNode;
	currDepth = -1;

	if (config.HEURST_FLAG == 1)
	{
		struct BnCnodeType *hrsticNode = NULL;  // heuristic root node
		hrsticNode = copyNode(rootNode, 0.1);
		hrsticNode->ishrstic = 1;
		rootNode->nextnode = hrsticNode;
		hrsticNode->prevnode = rootNode;
		activeNode = hrsticNode;
		currentNode = hrsticNode;
	}

	while (activeNode != NULL)
	{
		if (branchNode(stoc, prob, cell, currentNode, &activeNode))
			errMsg("BnC", "branchbound", "branching failed", 0);

		if (activeNode == NULL) break;
		if (activeNode->prevnode == NULL && nodecnt > 1) break;

		if (cell->k == config.MAX_ITER) break;
		if (nodecnt > config.MAX_NODES) break;

		currentNode = activeNode;
		nodecnt++;

		/*
		Loop for revisiting the fractional nodes to put them back to the queue if the estimate gets better
		*/
#if defined(useDNODE)
		if (dnodes > 0)
		{
			for (int cnt = 0; cnt < dnodes; cnt++)
			{
				double est = vXvSparse(nodearr[cnt]->vars, prob[0]->dBar) + maxCutHeight(cell->cuts, cell->sampleSize, nodearr[cnt]->vars, prob[0]->num->cols, prob[0]->lb);
				if (est < cell->incumbEst && est > meanVal)
				{
					nodearr[cnt]->isActive = true;
					for (int n = cnt; n < dnodes - 1; n++)
					{
						nodearr[n] = nodearr[n + 1];
					}
					nodearr[dnodes] = NULL;
					dnodes -= 1;
				}
			}
		}
#endif // defined(useDNODE)



		/*
		Loop for revisiting the integer feasible nodes to update the current best
		*/
#if defined(useINODE)
		if (inodes > 0)
		{
			for (int cnt = 0; cnt < inodes; cnt++)
			{
				double est = vXvSparse(inodearr[cnt]->vars, prob[0]->dBar) + maxCutHeight(cell->cuts, cell->sampleSize, inodearr[cnt]->vars, prob[0]->num->cols, prob[0]->lb);
				if (est < GlobeUB && est > meanVal)
				{
					GlobeUB = est;
					bestNode = inodearr[cnt];
				}
			}
		}
#endif // defined(useINODE)


	}

	//Replace the best node to the incumbent for the out of sample testing
	copyVector(bestNode->vars, cell->incumbX, bestNode->edInt, true);
	cell->incumbEst = bestNode->objval;

#if defined(printBest)
	printLine();
	printLine();
	printf("Best node key:  %d - depth: %d \n", bestNode->key, bestNode->depth);
	printf("Best Int feasible Solution: \n");
	printVector(bestNode->vars, bestNode->edInt, NULL);
	printf("\nNumber of disjunctions: %d \n", sumintVec(bestNode->disjncs, bestNode->edInt));
	printLine();
	printLine();
#endif // defined(printBest)

	//freeOneProblem(original);
	freeNodes(rootNode);
	return 0;

TERMINATE:
	//freeOneProblem(original);
	freeNodes(rootNode);
	return 1;
}

/* Return the previous active node */
struct BnCnodeType *nextNode(struct BnCnodeType *node)
{

	struct BnCnodeType *temp = NULL;
	temp = node;
	bool contin = true;

	if (node->isActive == false)
	{
		while (contin)
		{
			temp = temp->prevnode;

			if (temp->key == 0)
			{
				contin = false;
				break;
			}

			if (temp->isActive == true)
			{
				contin = false;
				break;
			}

		}
	}
	else
	{
		temp = node;

	}


	return temp;
}

/* get the index of the first leaf */
int getfirstLeaf(int depth)
{
	if (depth == 0) return 0;
	int out = 0;

	if (depth > 1)
	{
		for (int i = 1; i < depth; i++)
			out += pow(2, i);
	}


	return out + 1;
}

/* get the index of the node */
int getnodeIdx(int depth, int key, int isleft)
{
	int out = 0;

	out += getfirstLeaf(depth);
	out += (key - getfirstLeaf(depth - 1)) * 2;

	if (isleft == 0) out += 1;

	return out;
}


/* this function returns the fraction of lamda values are used from the earlier iterations or
from the tight lambdas of the parent and new lamdas discovered in the current node */
void fracLamda(cellType *cell, struct BnCnodeType *node)
{
	int i, j, k, notnewLambda, totLambda;

	node->Lambdasize = cell->lambda->cnt;
	notnewLambda = 0;
	totLambda = 0;

	if (cell->basis->basisEval == 1)
	{
		for (j = 0; j < cell->cuts->vals[cell->iCutIdx]->numSamples; j++)
		{
			if (cell->cuts->vals[cell->iCutIdx]->iStar[j] != -1)
			{
				node->IncumbiStar[totLambda] = cell->cuts->vals[cell->iCutIdx]->iStar[j];
				totLambda += 1;
				notnewLambda += 1;
			}
		}
	}


	node->tightPi = totLambda;
	node->fracPi = ((double)totLambda) / ((double)node->Lambdasize);


}

/* given a node problem this subroutine branch on the selected variable
it returns 0 if no branching is needed (solution of the node problem is integer)
or 1 when the node problem is fractional */
int branchNode(stocType *stoc, probType **prob, cellType *cell, struct BnCnodeType *node, struct BnCnodeType **activeNode)
{

	//Adjust the basis init if the pi eval flag is on
	if (node->key > 1 && config.Pi_EVAL_FLAG == 1)
	{
		cell->basis->init = cell->basis->cnt;
		for (int b = 0; b < maxcut*config.MAX_ITER; b++)
		{
			if (b < node->partightPi)
			{
				cell->basis->iStar = node->ParIncumbiStar[b];
			}
			else
			{
				cell->basis->iStar[b] = -1;
			}
		}
	}

	// Solve the node problem and obtain the solutions
	if (solveNode(stoc,prob,cell, node, original->name)) {

#if defined(printBranch)
		if (node->key > 0)
		{
			printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%-12s\n", node->key, node->parentkey, node->depth, cell->k, oneNorm(node->vars, node->edInt),0.0, 0.0, "False", "False", "NaN");
		}
		else
		{
			printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%-12s\n", node->key, 0, 0, cell->k, oneNorm(node->vars, node->edInt), 0.0, 0.0, "False", "False", "NaN");
		}
#endif // defined(printBranch)

		struct BnCnodeType * nodenew = NULL;
		if (node->prevnode != NULL)
		{
			nodenew = node;
			node = nodenew->prevnode;
			node->nextnode = NULL;
			*activeNode = nextNode(node);
		}
		else {
			errMsg("beanchBnC", "branchNode", "failed to solve the root node", 0);
			return 1;
		}

		freeNode(nodenew);
		return 0;
	}

	node->LB = node->objval;
	node->numSamp = cell->k;
	fracLamda(cell, node);

	// Check if the obtained solution from the solveNode is integer 
	if (isInteger(node->vars, node->edInt, 0, node->edInt + 1, config.TOLERANCE))
	{
#if defined(printBranch)
		if (node->key > 0)
		{
			printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%-10.3f\n", node->key, node->parentkey, node->depth, cell->k, oneNorm(node->vars, node->edInt), node->LB, GlobeUB, "True", "True",node->fracPi);
		}
		else
		{
			printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%-12s\n", node->key, 0, 0, cell->k, oneNorm(node->vars, node->edInt), node->LB, GlobeUB, "True", "True", "NaN");
		}
#endif // defined(printBranch)		
		
		if (node->prevnode == NULL) return 0;
		node->isActive = false;
		if (node->isSPopt == true && node->objval < GlobeUB )
		{
			GlobeUB = node->objval;
			bestNode = node;
		}
		else
		{
#if defined(useINODE)
			if (inodes < maxdnodes) inodearr[inodes++] = node;
#endif // defined(useINODE)
		}
		node->UB = node->objval;
		if (node->prevnode->depth == 0) *activeNode = NULL; else *activeNode = nextNode(node);
		currDepth = (*activeNode)->depth;

		return 0;
	}

#if defined(printBranch)
	if (node->key > 0)
	{
		printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%-10.3f\n", node->key, node->parentkey, node->depth, cell->k, oneNorm(node->vars, node->edInt), node->LB, GlobeUB, "True", "False",node->fracPi);
	}
	else
	{
		printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%-12s\n", node->key, 0, 0, cell->k, oneNorm(node->vars, node->edInt), node->LB, GlobeUB, "True", "False","NaN");
	}
#endif // defined(printBranch)


	node->UB = GlobeUB;

	// Compare the obj val with the Global LB
	if (node->objval > GlobeUB)
	{
		node->isActive = false;
#if defined(useDNODE)
		if (dnodes < maxdnodes) nodearr[dnodes++] = node;
#endif // defined(useINODE)
		if (node->prevnode->depth == 0) *activeNode = NULL; else *activeNode = nextNode(node);
		if (node->prevnode->depth == 0) currDepth = 0; else currDepth = (*activeNode)->depth;
		return 0;
	}

	// Incement depth
	currDepth = (*activeNode)->depth;

	//Select the variable index for branching 
	int vaIdx = branchVar(node, config.VAR_STR);

#if  defined(depthtest)
	if (node->depth == node->edInt)
	{
		printf("Solution at this depth: \n");
		for (int v = 0; v < node->edInt; v++)
		{
			printf("var%d:  val:  %0.04fl - LB: %0.04fl - UB: %0.04fl - distjnct: %d\n"
				, v, node->vars[v + 1], node->disjncsVal[v][0], node->disjncsVal[v][1], node->disjncs[v]);
		}
	}
#endif //  defined(depthtest)


	// Branch on the node
	node->isActive = false;
	if (node->depth < node->edInt)
	{
#if defined(BNC_CHECK)
		printf("\nbefore branching - varIdx %d:\n",vaIdx+1);
		printVector(node->vars, node->numVar, NULL);
#endif // defined(BNC_CHECK)
		node->right = newNode(getnodeIdx(node->depth + 1, node->key, 0), node, node->vars[vaIdx + 1], vaIdx, false);
		node->left = newNode(getnodeIdx(node->depth + 1, node->key, 1), node, node->vars[vaIdx + 1], vaIdx, true);
		node->right->prevnode = node;
		node->nextnode = node->right;
		node->right->nextnode = node->left;
		node->left->prevnode = node->right;
		*activeNode = node->left;
	}
	else
	{
		*activeNode = nextNode(node);
		currDepth = (*activeNode)->depth;
	}

	return 0;
}

