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

#undef writeprob

int currKey;
int currDepth;
double GlobeUB;             // Global upper bound
oneProblem      *original;  // Info of the original problem 
struct BnCnodeType *bestNode; // best node that is found so far
struct BnCnodeType **nodearr; // array of deactivated leaf nodes
int dnodes;                   /* number of deactivated nodes */
#define printBest
#define testBnC
#define printBranch
#undef depthtest
#undef printSol

int sumintVec(iVector a, int len)
{
	int n;
	int output = 0.0;

	for (n = 0; n < len; n++)
		output += a[n];

	return output;
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
	temp->left = temp->right = NULL;
	temp->nextnode = temp->prevnode = NULL;
	temp->numVar = numVar;
	temp->LB = LB;
	temp->UB = UB;
	temp->isActive = true;
	temp->isfathomed = false;
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
	if (fabs(fracVal - round(fracVal)) < 0.00001)
		fracVal = fracVal - 2 * config.TOLERANCE;
	struct BnCnodeType *temp = (struct BnCnodeType *)malloc(sizeof(struct BnCnodeType));
	temp->key = key;
	temp->parentkey = parent->key;
	temp->depth = parent->depth + 1;
	temp->isleft = isleft;
	temp->fracVal = fracVal;
	temp->varId = varId;
	temp->left = temp->right = NULL;
	temp->nextnode = temp->prevnode = NULL;
	temp->numVar = parent->numVar;
	temp->stInt = parent->stInt;
	temp->edInt = parent->edInt;
	temp->LB = parent->LB;
	temp->UB = parent->UB;
	temp->isActive = true;
	temp->isfathomed = false;
	if (!(temp->disjncs = (iVector)arr_alloc(temp->numVar, int)))
		errMsg("allocation", "newNode", "temp->disjncs", 0);
	if (!(temp->disjncsVal = (dVector *)arr_alloc(temp->numVar, dVector)))
		errMsg("allocation", "newNode", "temp->disjncs", 0);
	for (i = 0; i < parent->numVar; i++)
	{
		if (!(temp->disjncsVal[i] = (dVector)arr_alloc(2, double)))
			errMsg("allocation", "newNode", "temp->disjncs", 0);
		if (i == varId)
		{
			if (temp->isleft)
			{
				temp->disjncs[i] = 1;
				temp->disjncsVal[i][0] = parent->disjncsVal[i][0];
				temp->disjncsVal[i][1] = floor(fracVal);
			}
			else
			{
				temp->disjncs[i] = 1;
				temp->disjncsVal[i][0] = ceil(fracVal);
				temp->disjncsVal[i][1] = parent->disjncsVal[i][1];
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

	if (node->isleft)
	{
		node->vars[node->varId+1] = disjncsVal[node->varId][1];
	}
	else {
		node->vars[node->varId+1] = disjncsVal[node->varId][0];
	}

	
	if (changeQPbds(cell->master->lp, numCols, lbounds, ubounds, node->vars, 0)) {
		errMsg("algorithm", "algoIntSD", "failed to change the bounds to convert the problem into QP", 0);
		return 1;
	}

	mem_free(lbounds); 
	mem_free(ubounds);

	return 0;
}


 // Subroutine for Solving the node problem given the lp pointer and 
 // node information
double solveNode(stocType *stoc, probType **prob, cellType *cell, struct BnCnodeType *node, cString pname)
{
	int 	status;
	cell->ki = 0;

	if (node->prevnode != NULL)
	{
		if (addBnCDisjnct(cell, node->disjncsVal, node->edInt, node))
			errMsg("addDisjnct", "solveNode", "adding disjunctions are failed", 0);

		//Initializing candidX, candidEst, incumbEst and IncumbX
		copyVector(node->vars, cell->incumbX, node->edInt, true);
		copyVector(cell->incumbX, cell->candidX, node->edInt, true);
		cell->candidEst = vXvSparse(cell->candidX, prob[0]->dBar) + maxCutHeight(cell->cuts, cell->sampleSize, cell->candidX, prob[0]->num->cols, prob[0]->lb);
		cell->incumbEst = cell->candidEst;
		node->objval = cell->incumbEst;

		// Setup the QP master problem 
		if (changeQPproximal(cell->master->lp, node->edInt, cell->quadScalar)) {
			errMsg("algorithm", "algoIntSD", "failed to change the proximal term", 0);
			return 1;
		}

		if (changeQPrhs(prob[0], cell, node->vars)) {
			errMsg("algorithm", "algoIntSD", "failed to change the right-hand side to convert the problem into QP", 0);
			return 1;
		}

		if (cell->normDk >= config.R3 * cell->normDk_1) {
			cell->quadScalar *= config.R2 * config.R3 * cell->normDk_1 / cell->normDk;
			cell->quadScalar = minimum(config.MAX_QUAD_SCALAR, cell->quadScalar);
			cell->quadScalar = maximum(config.MIN_QUAD_SCALAR, cell->quadScalar);
		}

		/* update the candidate cut as the new incumbent cut */
		cell->iCutUpdt = cell->k;
		cell->incumbChg = true;

		/* keep the two norm of solution*/
		cell->normDk_1 = cell->normDk;
		/* Since incumbent solution is now replaced by a candidate, we assume it is feasible now */
		cell->infeasIncumb = false;
		/* gamma needs to be reset to 0 since there's no difference between candidate and incumbent*/
		cell->gamma = 0.0;
	}

#if defined(printSol)
	printVector(node->vars, node->edInt, NULL);
#endif // defined(printSol)
	
	if (node->prevnode == NULL || (node->objval < GlobeUB && !isInteger(node->vars, node->edInt, 0, node->edInt + 1, config.TOLERANCE)))
	{
		/* Use two-stage stochastic decomposition algorithm to solve the problem */
		if (solveCell(stoc, prob, cell)) {
			return 1;
		}

		copyVector(cell->incumbX, node->vars, node->edInt, true);
		node->objval = cell->incumbEst;
		cell->optFlag = false;
	}
	else
	{
		printf("-SD ends (dnodes=%d)\n", dnodes + 1);

		return 0;
	}

#if defined(printSol)
	printVector(node->vars, node->edInt, NULL);
#endif // defined(printSol)

	printf("-SD ends (dnodes=%d)\n",dnodes+1);

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

}


/* branch and bound algorithm */
int branchbound(stocType *stoc, probType **prob, cellType *cell, double LB, double UB)
{

	int i, j;

	if (!(nodearr = (struct BnCnodeType **)arr_alloc(maxdnodes, struct BnCnodeType *)))
		errMsg("allocation", "branchbound", "nodearr", 0);
	dnodes = -1;

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
	printf("%-10s%-10s%-10s%-10s%-10s%-12s%-12s%-12s%-12s%\n", "node id", "parent", "depth", "k", "\|x\|","fval", "UB", "feasible", "integer");
	printLongLine();
#endif // defined(printBranch)

	int totdepth = original->mac;
	double totnodes = pow(2, totdepth);

	rootNode = newrootNode(original->mac, LB, UB, original);
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
	activeNode = rootNode;
	currentNode = rootNode;
	currDepth = -1;

	while (activeNode != NULL)
	{
		if (branchNode(stoc, prob, cell, currentNode, &activeNode))
			errMsg("BnC", "branchbound", "branching failed", 0);

		if (activeNode->prevnode == NULL) break;

		if (cell->k == config.MAX_ITER) break;

		currentNode = activeNode;

		for (int cnt = 0; cnt < dnodes; cnt++)
		{
			double est = vXvSparse(nodearr[cnt]->vars, prob[0]->dBar) + maxCutHeight(cell->cuts, cell->sampleSize, cell->candidX, prob[0]->num->cols, prob[0]->lb);
			if (est < cell->incumbEst)
			{
				nodearr[cnt]->isActive = true;
				for (int n = cnt; n < dnodes-1; n++)
				{
					nodearr[n] = nodearr[n + 1];
				}
				nodearr[dnodes] = NULL;
				dnodes -= 1;
			}
		}
	}

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

	freeOneProblem(original);
	freeNodes(rootNode);
	return 0;

TERMINATE:
	freeOneProblem(original);
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

/* given a node problem this subroutine branch on the selected variable
it returns 0 if no branching is needed (solution of the node problem is integer)
or 1 when the node problem is fractional */
int branchNode(stocType *stoc, probType **prob, cellType *cell, struct BnCnodeType *node, struct BnCnodeType **activeNode)
{


	// Solve the node problem and obtain the solutions
	if (solveNode(stoc,prob,cell, node, original->name)) {

#if defined(printBranch)
		if (node->key > 0)
		{
			printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%\n", node->key, node->parentkey, node->depth, cell->k, oneNorm(node->vars, node->edInt),0.0, 0.0, "False", "False");
		}
		else
		{
			printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%\n", node->key, 0, 0, cell->k, oneNorm(node->vars, node->edInt), 0.0, 0.0, "False", "False");
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


	// Check if the obtained solution from the solveNode is integer 
	if (isInteger(node->vars, node->edInt, 0, node->edInt + 1, config.TOLERANCE))
	{
#if defined(printBranch)
		if (node->key > 0)
		{
			printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%\n", node->key, node->parentkey, node->depth, cell->k, oneNorm(node->vars, node->edInt), node->LB, GlobeUB, "True", "True");
		}
		else
		{
			printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%\n", node->key, 0, 0, cell->k, oneNorm(node->vars, node->edInt), node->LB, GlobeUB, "True", "True");
		}
#endif // defined(printBranch)		
		
		if (node->prevnode == NULL) return 0;
		node->isActive = false;
		if (node->objval < GlobeUB)
		{
			GlobeUB = node->objval;
			bestNode = node;
		}
		else
		{
			if (dnodes < maxdnodes) nodearr[dnodes++] = node;
		}
		node->UB = node->objval;
		if (node->prevnode->depth == 0) *activeNode = NULL; else *activeNode = nextNode(node);
		currDepth = (*activeNode)->depth;

		return 0;
	}

#if defined(printBranch)
	if (node->key > 0)
	{
		printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%\n", node->key, node->parentkey, node->depth, cell->k, oneNorm(node->vars, node->edInt), node->LB, GlobeUB, "True", "False");
	}
	else
	{
		printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%\n", node->key, 0, 0, cell->k, oneNorm(node->vars, node->edInt), node->LB, GlobeUB, "True", "False");
	}
#endif // defined(printBranch)


	node->UB = GlobeUB;

	// Compare the obj val with the Global LB
	if (node->objval > GlobeUB)
	{
		node->isActive = false;
		if (dnodes < maxdnodes) nodearr[dnodes++] = node;
		if (node->prevnode->depth == 0) *activeNode = NULL; else *activeNode = nextNode(node);
		currDepth = (*activeNode)->depth;
		return 0;
	}

	// Incement depth
	currDepth = (*activeNode)->depth;

	//Select the variable index for branching 
	int vaIdx = branchVar(node, 0);

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
