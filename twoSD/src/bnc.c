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

#include "bnc.h"


extern configType config;


/* branch and bound algorithm */
int branchbound(stocType *stoc, probType **prob, cellType *cell, double LB, double UB) {
	int i;
	int nodecnt = 0;
	int repeat  = 0;
	int check_node_num = 0;
	int maxcut = config.CUT_MULT * cell->master->mac + 3;
	cell->basis->incumPicnt = maxcut;
	cell->basis->basisEval = config.Pi_EVAL_FLAG;
	cell->maxiter_nodes = 0;

#if defined(useDNODE)
	if (!(nodearr = (struct BnCnodeType **)arr_alloc(maxdnodes, struct BnCnodeType *)))
		errMsg("allocation", "branchbound", "nodearr", 0);
	if (!(inodearr = (struct BnCnodeType **)arr_alloc(maxdnodes, struct BnCnodeType *)))
		errMsg("allocation", "branchbound", "inodearr", 0);
#endif // defined(useDNODE)

	dnodes = -1;
	inodes = -1;

	/* set of LB and UB */
	GlobeUB = INFINITY;
	GlobeLB = INFINITY;

	original = prob[0]->sp;
	rootNode = NULL;  // root node
	bestNode = NULL;
	struct BnCnodeType * activeNode = NULL;// active node in the tree
	struct BnCnodeType * currentNode = NULL;// current node that is investigated

#if defined(printBranch)
	printLine();
	printf("Starting the BnC branch and bound procedure...\n");
	printf("Github page: https://github.com/siavashtab \n");
	printf("Email: stabrizian@smu.edu\n");
	printLine();
	printf("\n\n");
	printLine();
	printf("%-10s%-10s%-10s%-10s%-10s%-12s%-12s%-12s%-12s%-12s\n", "node id", "parent", "depth", "k", "\\|x\\|","fval", "UB", "feasible", "integer","Lamfrac");
	printLine();
#endif // defined(printBranch)

	rootNode = newrootNode(prob[0]->num->cols+1, LB, UB, original);
	rootNode->parobjVal = prob[0]->lb;
	meanVal = prob[0]->lb;
	cell->incumbEst = meanVal;
	cell->candidEst = meanVal;
	cell->tot_nodes = 1;
	cell->depth = 0;
	cell->maxiter_nodes = 0;

	for (i = 0; i < prob[0]->num->cols; i++) {
		if (original->ctype[i] == 'B' || original->ctype[i] == 'I') {
			rootNode->stInt = i;
			break;
		}
	}
	rootNode->edInt = 1;
	for (i = rootNode->stInt + 1; i < prob[0]->num->cols; i++) {
		if (original->ctype[i] != 'B' && original->ctype[i] != 'I') {
			rootNode->edInt = i;
			break;
		}
		else {
			rootNode->edInt++;
		}
	}
	rootNode->numIntVar = rootNode->edInt - rootNode->stInt;

	/* Copy the incumbent solution (currently stores the mean value problem solution) to node->vars. */
	for (i = 0; i < prob[0]->num->cols+1; i++)
		rootNode->vars[i] = cell->incumbX[i];
	activeNode = rootNode;
	currentNode = rootNode;
	currDepth = -1;

	if (config.HEURST_FLAG == 1) {
		struct BnCnodeType *hrsticNode = NULL;  // heuristic root node
		hrsticNode = copyNode(rootNode, 0.1);
		hrsticNode->ishrstic = 1;
		rootNode->nextnode = hrsticNode;
		hrsticNode->prevnode = rootNode;
		activeNode = hrsticNode;
		currentNode = hrsticNode;
	}

	/* Loop though all the active nodes (leaf nodes) to update their incumbent estimates */
	while (activeNode != NULL) {
//	while ( nodecnt < 2 ) {

		if (branchNode(stoc, prob, cell, currentNode, &activeNode))
			errMsg("BnC", "branchbound", "branching failed", 0);


		//Revising the paused nodes
#if defined(useDNODE)
		/* Loop for revisiting the fractional nodes to put them back to the queue if the estimate gets better */
		for (int cnt = 0; cnt < dnodes; cnt++) {
			if (cell->cutsPool[nodearr[cnt]->poolID]) {
				for (int cnt2 = 0; cnt2 < cell->cutsPool[nodearr[cnt]->poolID]->cnt; cnt2++)
					revisitNode(prob[1]->num, prob[1]->coord, cell->basis, cell->sigma, cell->delta, cell->omega, cell->sample,
						nodearr[cnt]->vars, cell->sampleSize, &cell->dualStableFlag, cell->pi_ratio, cell->k, cell->lb, cell->cutsPool[nodearr[cnt]->poolID]->vals[cnt2]);
				double est = vXvSparse(nodearr[cnt]->vars, prob[0]->dBar) +
					maxCutHeight(cell->cutsPool[nodearr[cnt]->poolID], cell->sampleSize, nodearr[cnt]->vars, prob[0]->num->cols, prob[0]->lb);
				if (est < cell->incumbEst && est > meanVal) {
					nodearr[cnt]->isActive = true;
					for (int n = cnt; n < dnodes - 1; n++) {
						nodearr[n] = nodearr[n + 1];
					}
					nodearr[dnodes] = NULL;
					dnodes -= 1;
				}
			}
		}
#endif // defined(useDNODE)

#if defined(useINODE)
		/* Loop for revisiting the integer feasible nodes to update the current best */
		for (int cnt = 0; cnt < inodes; cnt++) {
			if(cell->cutsPool[inodearr[cnt]->poolID]){
				for (int cnt2 = 0; cnt2 < cell->cutsPool[inodearr[cnt]->poolID]->cnt; cnt2++)
					revisitNode(prob[1]->num, prob[1]->coord, cell->basis, cell->sigma, cell->delta, cell->omega, cell->sample,
						inodearr[cnt]->vars, cell->sampleSize, &cell->dualStableFlag, cell->pi_ratio, cell->k, cell->lb, cell->cutsPool[inodearr[cnt]->poolID]->vals[cnt2]);
				double est = vXvSparse(inodearr[cnt]->vars, prob[0]->dBar)
					+ maxCutHeight(cell->cutsPool[inodearr[cnt]->poolID], cell->sampleSize, inodearr[cnt]->vars, prob[0]->num->cols, prob[0]->lb);
				if (est > inodearr[cnt]->parobjVal && est < GlobeUB && est > meanVal) {
					GlobeUB = est;
					bestNode = inodearr[cnt];
				}
			}
		}
#endif // defined(useINODE)

		//Checking the termination criteria
		if (activeNode == NULL) break;
		if (activeNode->prevnode == NULL && nodecnt > 3) break;
		if (activeNode->key == 0 && nodecnt > 3) break;
		if(activeNode->prevnode) if (activeNode->prevnode->key == 0 && nodecnt > 3) break;
		if (cell->k == config.MAX_ITER) break;
		if (nodecnt > config.MAX_NODES) break;
		if (repeat > 2) break;

		currentNode = activeNode;
		nodecnt = cell->tot_nodes;

		if (nodecnt > check_node_num) {
			check_node_num = nodecnt;
		}
		else {
			repeat++;
		}

	}

	if(inodes > 0) cell->int_nodes = inodes;
	if(dnodes > 0) cell->d_nodes = dnodes;

	if ( bestNode != NULL ) {

		if (bestNode->objval < GlobeLB) {
			// Replace the best node to the incumbent for the out-of-sample testing
			copyVector(bestNode->vars, cell->incumbX, bestNode->numVar, true);
			cell->incumbEst = bestNode->objval;

			// Remove the tolerance from the best solution
			for (int v = 0; v < bestNode->edInt; v++)
				bestNode->vars[v] = round(bestNode->vars[v]);
		}
		else {
			goto USEGLOBALLB;
		}
	}
	else {

	USEGLOBALLB:

		// Remove the tolerance from the best solution
		for (int v = 0; v < prob[0]->num->cols+1; v++)
			cell->incumbX[v] = LBround[v];

		// Update the incumbent estimate 
		cell->incumbEst = GlobeLB;
	}

	/* cleaning the master problem after finishing bnc */
	if (cleanMaster(prob[0], cell)) {
		errMsg("BnC", "solveNode", "failed to clean the master problem after finishing bnc ", 0);
		return 1;
	}

#if defined(printBest)
	printLine();
	printLine();
	printf("Best node key:  %d - depth: %d \n", bestNode->key, bestNode->depth);
	printf("Best Int feasicopy->name = (cString) arr_alloc(NAMESIZE, char);ble Solution: \n");
	printVector(bestNode->vars, bestNode->edInt, NULL);
	printf("\nNumber of disjunctions: %d \n", sumintVec(bestNode->disjncs, bestNode->edInt));
	printLine();
	printLine();
#endif // defined(printBest)

	// Free node queue
	if (rootNode->nextnode && rootNode->nextnode->key > 0) {
		freeNodes(rootNode);
	}
	else {
		if(rootNode->disjncs) freeNode(rootNode);
	}
	if (LBround)  mem_free(LBround);
	mem_free(nodearr);
	mem_free(inodearr);

	return 0;
}//END branchBound()

/* given a node problem this subroutine branch on the selected variable
it returns 0 if no branching is needed (solution of the node problem is integer)
or 1 when the node problem is fractional */
int branchNode(stocType *stoc, probType **prob, cellType *cell, struct BnCnodeType *node, struct BnCnodeType **activeNode) {

	/* Identify the bases generated by the parent node. This is used to streamline the argmax operation to explore only
	the bases generated by the parent as all the bases may not useful. */
	if (node->key > 1 && config.Pi_EVAL_FLAG == 1) {
		cell->basis->init = cell->basis->cnt;
		for (int b = 0; b < maxcut*config.MAX_ITER; b++) {
			if (b < node->partightPi) {
				cell->basis->iStar[b] = node->ParIncumbiStar[b];
			}
			else {
				cell->basis->iStar[b] = -1;
			}
		}
	}

	// Solve the node problem and obtain the solutions
	int status = solveNode(stoc, prob, cell, node, original->name);
	int key = node->key;
	

	if ( !cell->masterFeasFlag ) {
#if defined(printBranch)
		if (node->key > 0) {
			printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%-12s\n",
					node->key, node->parentkey, node->depth, cell->k, oneNorm(node->vars, node->edInt),0.0, 0.0, "False", "False", "NaN");
		}
		else {
			printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%-12s\n",
					node->key, 0, 0, cell->k, oneNorm(node->vars, node->edInt), 0.0, 0.0, "False", "False", "NaN");
		}
#endif // defined(printBranch)

		if (node->prevnode != NULL) {
			/* Active a new node */
			*activeNode = nextNode(node);
		}
		else {
			errMsg("beanchBnC", "branchNode", "failed to solve the root node", 0);
			return 1;
		}

		freepartialNode(node);
		return 0;
	}
	else if (status == 1) {
		if (cell->ki == config.MAX_ITER_CLBK) {
			if (node->depth == node->edInt || isInteger(node->vars, node->numIntVar, node->stInt, node->edInt, config.INT_TOLERANCE))
				goto INTEGER;
		}
		if (node->prevnode != NULL) {
			*activeNode = nextNode(node);
			freepartialNode(node);
		}
		else {
			errMsg("beanchBnC", "branchNode", "failed to solve the root node", 0);
			return 1;
		}

		return 0;
	}

	node->LB = node->objval;
	node->numSamp = cell->k;

	// Condition 1. Check if the obtained solution from the solveNode is integer
	if (node->depth == node->edInt || isInteger(node->vars, node->numIntVar, node->stInt, node->edInt, config.INT_TOLERANCE)) {
		INTEGER:
#if defined(printBranch)
		if (node->key > 0) {
			printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%-10.2f\n",
					node->key, node->parentkey, node->depth, cell->k, oneNorm(node->vars, node->edInt), node->LB, GlobeUB, "True", "True",node->fracPi);
		}
		else {
			printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%-12s\n",
					node->key, 0, 0, cell->k, oneNorm(node->vars, node->edInt), node->LB, GlobeUB, "True", "True", "NaN");
		}
#endif // defined(printBranch)

		if (node->prevnode == NULL)
			return 0;
		node->isActive = false;

		if (node->isSPopt == true && node->objval < GlobeUB ) {
			GlobeUB = node->objval;
			bestNode = node;
		}
		else {
#if defined(useINODE)
			if (inodes < maxdnodes)
				inodearr[++inodes] = node;
#endif // defined(useINODE)
		}
		node->UB = node->objval;
		if (node->prevnode->depth == 0) {
			return 0;
		}
		else {
			*activeNode = nextNode(node);
		}
		currDepth = (*activeNode)->depth;

		return 0;
	}

#if defined(printBranch)
	if (node->key > 0) {
		printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%-10.2f\n", node->key, node->parentkey, node->depth, cell->k, oneNorm(node->vars, node->edInt), node->LB, GlobeUB, "True", "False",node->fracPi);
	}
	else {
		printf("%-10d%-10d%-10d%-10d%-10.1f%-14.2f%-14.2f%-12s%-12s%-12s\n", node->key, 0, 0, cell->k, oneNorm(node->vars, node->edInt), node->LB, GlobeUB, "True", "False","NaN");
	}
#endif // defined(printBranch)

	node->UB = GlobeUB;

	// Condition 2. The obtained solution is below the global lower bound, add the node to a pool.
	// Compare the obj val with the Global LB
	if ((node->objval > GlobeUB && node->key > 2) || (node->isSPopt == false && bestNode != NULL)) {
		if (node->prevnode) {
			node->isActive = false;
#if defined(useDNODE)
			if (dnodes < maxdnodes)
				nodearr[++dnodes] = node;
#endif // defined(useINODE)
			if (node->prevnode->key == 0) {
				*activeNode = NULL;
			}
			else {
				*activeNode = nextNode(node);
			}

			if (node->prevnode->depth == 0) currDepth = 0; else currDepth = (*activeNode)->depth;
			return 0;
		}

	}

	// Condition 3. Select the variable index for branching
	currDepth = (*activeNode)->depth;
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
	struct BnCnodeType *right, *left;

	right = newNode(getnodeIdx(node->depth + 1, node->key, 0), node, node->vars[vaIdx + 1], vaIdx, false);
	right->poolID = node->poolID;

	left  = newNode(getnodeIdx(node->depth + 1, node->key, 1), node, node->vars[vaIdx + 1], vaIdx, true);
	left->poolID = cell->numPools++;
	cell->cutsPool[left->poolID] = NULL;
	if (cell->cutsPool[left->parentPoolID] != NULL && cell->cutsPool[left->parentPoolID]->cnt > 0) {
		copyCuts(prob[0]->num, cell->cutsPool[left->parentPoolID], &cell->cutsPool[left->poolID]);
	}
	else {
		printf("parent node of %d cut pool is empty", left->key);
	}
	if (left->depth > cell->depth) cell->depth = left->depth;
	cell->tot_nodes += 2;


	if ( node->key == 0 ) {
		node->nextnode = right;
		right->nextnode = left;

		left->prevnode = right;
		right->prevnode = node;
	}
	else {

		if (node->nextnode) {
			node->nextnode->prevnode = left;
			left->nextnode = node->nextnode;
		}

		node->nextnode = right;
		right->prevnode = node;

		left->prevnode = right;
		right->nextnode = left;

		freepartialNode(node);
	}

	*activeNode = left;


	return 0;
}//END branchNode()

/* given a node problem which is solved decide which variable should be selected
 for branching  */
int branchVar(struct BnCnodeType *node, int strategy) {

	int idx, maxidx, minidx;
	double larger, smaller, avg, frac;
	
	if (strategy == 1)// strategy based on the larger fractional value
	{
		larger = -1;  /* maximum fractional value */
		maxidx = node->depth;

		for (int i = 0; i < node->numIntVar; i++)
		{
			if (node->disjncs[i] != 1 && larger < node->vars[i+1])
			{
				larger = node->vars[i+1];
				maxidx = i;
			}
		}

		return maxidx;
	}
	else if (strategy == 2)// strategy based on the smaller fractional value
	{
		smaller = INFINITY;  /* minimum fractional value */
		minidx = node->depth;

		for (int i = 0; i < node->numIntVar; i++)
		{
			if (node->disjncs[i] != 1 && smaller > node->vars[i+1])
			{
				smaller = node->vars[i+1];
				minidx = i;
			}
		}

		return minidx;
	}
	else if (strategy == 3)// strategy based on the value closer to 0.5
	{
		avg = INFINITY;  /* avg fractional value */
		idx = node->depth;

		for (int i = 0; i < node->numIntVar; i++)
		{
			if (node->disjncs[i] != 1 && abs(0.5 - node->vars[i+1]) < avg)
			{
				avg = abs(0.5 - node->vars[i+1]);
				idx = i;
			}
		}

		return idx;
	}

	idx = node->depth;

	for (int i = 0; i < node->numIntVar; i++)
	{
		
		frac = fabs(node->vars[i + 1] - round(node->vars[i + 1]));

		if (node->disjncs[i] != 1 && frac > config.INT_TOLERANCE)
		{
			idx = i;
			return i;
		}
	}

	return idx;
}//END branchVar()

// Subroutine for Solving the node problem given the lp pointer and node information
int solveNode(stocType *stoc, probType **prob, cellType *cell, struct BnCnodeType *node, cString pname) {
	cell->ki = 0;

#if defined(BNC_CHECK)
	printf("\ndisjuction: ");
	for (int n = 0; n < node->numIntVar; n++)
		printf("%d, ", node->disjncs[n]);
	printf("\n");
#endif // defined(BNC_CHECK)

	/* 1. Setup the node to be solved */
	if ( setupNode(prob[0], cell, node) ) {
		errMsg("BnC", "solveNode", "failed to setup node before solveNode", 0);
		return 1;
	}

#if defined(BNC_CHECK)
	printf("\ninit var: ");
	printVector(node->vars, node->numVar-1, NULL);
	printf("\ninit incumbX: ");
	printVector(cell->incumbX, node->numVar-1, NULL);
	printf("\ninit candidX: ");
	printVector(cell->candidX, node->numVar-1, NULL);
#endif // defined(BNC_CHECK)

	/* 2. Invoke the SD solver to solve the node */
	if (node->ishrstic ||
		(cell->k < config.MAX_ITER && (node->prevnode == NULL ||
		(node->objval < GlobeUB && sumDisjncs(node->disjncs, node->numIntVar)))) &&
		node->lbdsjncs < 1*node->numIntVar && node->ubdsjncs < 0.8*node->numIntVar) {

		/* solve the current master problem */
		int status;
		changeQPSolverType(ALG_CONCURRENT);
		if (solveProblem(cell->master->lp, cell->master->name, cell->master->type, &status)) {
			errMsg("BnC", "solveNode", "failed to solve the node using SD", 0);
			if (cleanInfNode(prob[0], cell, node)) {
				errMsg("BnC", "solveNode", "failed to clean the inf node after SD solve", 0);
				return 1;
			}
#if defined(BNC_INF_CHECK)
			char fname[NAMESIZE];
			sprintf(fname, "%s_n%d_k%d.lp", "master", node->key, cell->k);
			writeProblem(cell->master->lp, fname);
#endif // defined(BNC_INF_CHECK)
			return 1;
		}

		/* Use two-stage stochastic decomposition algorithm to solve the problem */
		bool isRoot = node->key == 0;
		if ( solveCell(stoc, prob, cell, isRoot) ) {
			errMsg("BnC", "solveNode", "failed to solve the node using SD", 0);
			if (cleanInfNode(prob[0], cell, node)) {
				errMsg("BnC", "solveNode", "failed to clean the inf node after SD solve", 0);
				return 1;
			}
#if defined(BNC_INF_CHECK)
			char fname[NAMESIZE];
			sprintf(fname, "%s_n%d_k%d.lp", "master",node->key, cell->k);
			writeProblem(cell->master->lp, fname);
#endif // defined(BNC_INF_CHECK)

			return 1;
		}
		node->isSPopt = true;


		// If the solution is fractional create a rounded solution and update the GlobalLB
		if (!isInteger(node->vars, node->numIntVar, node->stInt, node->edInt, config.INT_TOLERANCE)) {
			for (int v = 0; v < node->numVar; v++) {
				LBround[v] = round(cell->incumbX[v]);
			}

			double LBlocal = vXvSparse(LBround, prob[0]->dBar)
				+ maxCutHeight(cell->activeCuts, cell->sampleSize, LBround, prob[0]->num->cols, prob[0]->lb);
			
			if (LBlocal < GlobeLB) GlobeLB = LBlocal;
		}

#if defined(writemaster)
		char mname[NAMESIZE];
		sprintf(mname, "%s_k%d_n%d.lp", "masteraftersolve", cell->k, node->key);
		writeProblem(cell->master->lp, mname);
#endif // defined(writemaster)

		/* 2b. Clean the node before exit. */
		if ( cleanNode(prob[0], cell, node) ) {
			errMsg("BnC", "solveNode", "failed to clean the node after SD solve", 0);
			return 1;
		}

		if (cell->ki >= config.MAX_ITER_CLBK) {
			node->isSPopt = false;
			cell->maxiter_nodes++;
		}

	}
	else {
		truncate(node->vars, prob[0]->sp->bdl, prob[0]->sp->bdu, node->numVar);
		if (cell->incumbEst <= node->parobjVal || cell->incumbEst <= meanVal || cell->ki >= config.MAX_ITER_CLBK) {
			node->objval = meanVal - fabs(cell->incumbEst <= meanVal);
			node->isSPopt = false;
		}
		else {
			node->objval = cell->incumbEst;
		}
	}

	if (cell->ki > 600) printf("\n");

	printf("\nSD output: iters:%-4d - dnodes:%-3d - inodes:%-3d - sigma size:%-7d - lambda size:%-7d - omega size:%-7d\n",
			cell->ki, dnodes + 1, inodes +1, cell->sigma->cnt, cell->lambda->cnt, cell->omega->cnt);

	return 0;
}//END solveNode()

/* The subroutine sets up the B&B node by updating the cell structure with information necessary to solve the node SP. */
int setupNode(probType *prob, cellType *cell, struct BnCnodeType *node) {

	if ( node->key != 0) {

		if (cell->cutsPool[node->poolID] != NULL && cell->cutsPool[node->poolID]->cnt > 0) {
			/* 1a. Copy active cuts corresponding to the parent node from the cuts pool */
			copyCuts(prob->num, cell->cutsPool[node->poolID], &cell->activeCuts);

			/* 1b. Add the active cuts to the master problem */
			for (int n = 0; n < cell->activeCuts->cnt; n++) {
				if (addCut2Master(cell->master, cell->activeCuts->vals[n], cell->incumbX, prob->num->cols, false)) {
					errMsg("algoIntSD", "setupNode", "failed to add the new cut to master problem", 0);
					return 1;
				}
				cell->activeCuts->vals[n]->rowNum = prob->num->rows + n;
			}
		}

		/* 2. Retrieve the incumbent from the parent node and perform updates. */
		/* 2a. Add the B&B conditions */
		if (addBnCDisjnct(cell, node->disjncsVal, node,prob[0].sp->bdl,prob[0].sp->bdu)) {
			errMsg("addDisjnct", "solveNode", "adding disjunctions are failed", 0);
			return 1;
		}

		copyVector(node->vars, cell->incumbX, node->numVar, true);
		copyVector(cell->incumbX, cell->candidX, node->numVar, true);

		cell->candidEst = vXvSparse(cell->candidX, prob->dBar)
									+ maxCutHeight(cell->activeCuts, cell->sampleSize, cell->candidX, prob->num->cols, prob->lb);
		cell->incumbEst = cell->candidEst;
		node->objval = cell->incumbEst;
		cell->spFeasFlag = true;
	}

	/* 2b. Setup the quadratic master using the new incumbent */
	if (changeQPrhs(prob, cell, cell->incumbX)) {
		errMsg("algoIntSD", "setupNode", "failed to change the right-hand side to convert the problem into QP", 0);
		return 1;
	}

	/* 2c. Update the proximal parameter */
	cell->quadScalar = config.MIN_QUAD_SCALAR;
	if (changeQPproximal(cell->master->lp, prob->num->cols, cell->quadScalar)) {
		errMsg("algoIntSD", "setupNode", "failed to change the proximal term", 0);
		return 1;
	}

	/* 3. nodify the eta column based on the most recent number of observations */
	if (changeEtaCol(cell->master->lp, prob->num->rows, prob->num->cols, 1.0, cell->activeCuts)) {
		errMsg("algoIntSD", "setupNode", "failed to change the eta column coefficients", 0);
		return 1;
	}

	node->isSPopt = true;

#if defined(writemaster)
	char mname[NAMESIZE];
	sprintf(mname, "%s_k%d_n%d.lp", "masteraftersetup", cell->k, node->key);
	writeProblem(cell->master->lp, mname);
#endif // defined(writemaster)

	return 0;
}//END setupNode()

int cleanNode(probType *prob, cellType *cell, struct BnCnodeType *node) {

	/* 1. Copy the incumbent solution and estimate to the ndde structure */
	truncate(cell->incumbX, prob->sp->bdl, prob->sp->bdu, node->numIntVar);
	copyVector(cell->incumbX, node->vars, node->numVar, true);
	if ((cell->incumbEst <= 0.4*node->parobjVal && cell->incumbEst>=0) || 
		(abs(cell->incumbEst) >= abs(node->parobjVal) && cell->incumbEst <= 0) || 
		cell->incumbEst <= meanVal) {
		node->objval = meanVal - fabs(cell->incumbEst <= meanVal);
		node->isSPopt = false;
	}
	else {
		node->objval = cell->incumbEst;
	}
	cell->optFlag = false;

	// Update the info about the fraction of bases used in the incumbent cut for this node
	fracLamda(cell, node);

#if defined(BNC_CHECK)
	printf("\nafter SD var: ");
	printVector(node->vars, node->numVar-1, NULL);
#endif // defined(BNC_CHECK)

	/* 2. Copy the active cuts to the cutsPool */
	copyCutstoNodePool(prob->num, cell->activeCuts, cell->cutsPool[node->poolID], cell->piM);

	/* 3. Clean the master by removing all the inactive cuts */
	if (prob->num->rows < cell->master->mar) {
		for (int cnt = cell->master->mar - 1; cnt >= prob->num->rows; cnt--)
			if (removeRow(cell->master->lp, cnt, cnt)) {
				printf("row Num %d - tot rows %d - orig rows %d", cnt, cell->master->mar, prob->num->rows);
				errMsg("solver", "cleanNode", "failed to remove a row from master problem", 0);
				return 1;
			}
		cell->master->mar -= cell->activeCuts->cnt;
	}

	/* 4. Remove the all the cuts from the activeCuts structure */
	freeCutsType(cell->activeCuts, true);

#if defined(writemaster)
	char mname[NAMESIZE];
	sprintf(mname, "%s_k%d_n%d.lp", "masterafterclean", cell->k, node->key);
	writeProblem(cell->master->lp, mname);
#endif // defined(writemaster)
	return 0;
}//END cleanNode()

int cleanInfNode(probType *prob, cellType *cell, struct BnCnodeType *node) {

	int numCols = node->numVar - 1; /* eta column excluded */

	/* 0. Update the flags */
	node->isSPopt = false;
	node->isActive = false;
	node->isSPopt = false;
	

	/* 1. Clean the master by removing all the inactive cuts */
	if (prob->num->rows < cell->master->mar) {
		for (int cnt = cell->master->mar - 1; cnt >= prob->num->rows; cnt--)
			if (removeRow(cell->master->lp, cnt, cnt)) {
				printf("row Num %d - tot rows %d - orig rows %d", cnt, cell->master->mar, prob->num->rows);
				errMsg("solver", "cleanNode", "failed to remove a row from master problem", 0);
				return 1;
			}
		cell->master->mar -= cell->activeCuts->cnt;
	}

#if defined(BNC_CHECK)
	printf("\nafter SD var: ");
	printVector(node->vars, node->numVar, NULL);
#endif // defined(BNC_CHECK)


#if defined(writemaster)
	char mname[NAMESIZE];
	sprintf(mname, "%s_k%d_n%d.lp", "masterafterclean", cell->k, node->key);
	writeProblem(cell->master->lp, mname);
#endif // defined(writemaster)

	return 0;
}//END cleanInfNode()

/* Cleaning the master problem prob[0] after getting done with the branch and bound */
int cleanMaster(probType *prob, cellType *cell) {

	int 	cnt;
	int numCols = prob->num->cols; /* eta column excluded */

	/* Update and the bounds of variables */
	for (cnt = 0; cnt < numCols; cnt++) {
		cell->master->bdl[cnt] = prob->sp->bdl[cnt];
		cell->master->bdu[cnt] = prob->sp->bdu[cnt];
		//printf("\nv:%d  -  lb:%0.4f  -  ub:%0.4f", cnt, lbounds[cnt], ubounds[cnt]);
	}

	return 0;
}//END cleanMaster()

struct BnCnodeType *newrootNode(int numVar, double LB, double UB, oneProblem * orig) {
	int i;

	struct BnCnodeType *temp = (struct BnCnodeType *)malloc(sizeof(struct BnCnodeType));
	temp->key = 0;
	temp->parentkey = -1;
	temp->depth = 0;
	temp->isleft = false;
	temp->fracVal = -1;
	temp->varId = -1;
	temp->ishrstic = false;
	temp->nextnode = temp->prevnode = NULL;
	temp->numVar = numVar;
	temp->numIntVar = numVar - 1;
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
	temp->parentPoolID = -1;
	temp->poolID = 0;
	temp->lbdsjncs = 0;
	temp->ubdsjncs = 0;

	if (!(temp->disjncs = (iVector)arr_alloc(temp->numVar-1, int)))
		errMsg("allocation", "newNode", "temp->disjncs", 0);
	if (!(temp->disjncsVal = (dVector *)arr_alloc(temp->numVar-1, dVector)))
		errMsg("allocation", "newNode", "temp->disjncs", 0);
	for (i = 0; i < numVar-1; i++)
	{
		if (!(temp->disjncsVal[i] = (dVector)arr_alloc(2, double)))
			errMsg("allocation", "newNode", "temp->disjncs", 0);
	}
	if (!(temp->vars = (dVector)arr_alloc(temp->numVar, double)))
		errMsg("allocation", "newNode", "temp->vars", 0);
	if (!(LBround = (dVector)arr_alloc(temp->numVar, double)))
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
	else
	{
		temp->IncumbiStar = NULL;
		temp->ParIncumbiStar = NULL;
	}
	for (int v = 0; v < temp->numVar; v++)
		temp->vars[v] = 0.0;
	for (i = 0; i < temp->numVar-1; i++)
	{
		temp->disjncsVal[i][0] = orig->bdl[i];
		temp->disjncsVal[i][1] = orig->bdu[i];
		temp->disjncs[i] = 0;
	}

	return temp;

}//End newrootNode()

struct BnCnodeType *newNode(int key, struct BnCnodeType * parent, double fracVal, int varId, bool isleft) {
	int i;
	if (fabs(fracVal - round(fracVal)) == 0) {
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
	temp->nextnode = temp->prevnode = NULL;
	temp->numVar = parent->numVar;       /* numVar has oneNorm in the first element - sometimes it can be also seen as the 
										    total number of columns including eta */
	temp->numIntVar = parent->numIntVar;
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
	temp->parentPoolID = parent->poolID;
	temp->lbdsjncs = parent->lbdsjncs;
	temp->ubdsjncs = parent->ubdsjncs;

	if (!(temp->disjncs = (iVector)arr_alloc(temp->numIntVar, int)))
		errMsg("allocation", "newNode", "temp->disjncs", 0);
	if (!(temp->disjncsVal = (dVector *)arr_alloc(temp->numIntVar, dVector)))
		errMsg("allocation", "newNode", "temp->disjncs", 0);
	if (config.Pi_EVAL_FLAG == 1)
	{
		if (!(temp->IncumbiStar = (iVector)arr_alloc(config.MAX_ITER, int)))
			errMsg("allocation", "newNode", "temp->vars", 0);
		if (!(temp->ParIncumbiStar = (iVector)arr_alloc(config.MAX_ITER, int)))
			errMsg("allocation", "newNode", "temp->vars", 0);
		for (i = 0; i < config.MAX_ITER; i++)
		{
			temp->IncumbiStar[i] = -1;
			temp->ParIncumbiStar[i] = parent->IncumbiStar[i];
		}
	}
	else
	{
		temp->IncumbiStar = NULL;
		temp->ParIncumbiStar = NULL;
	}


	for (i = 0; i < parent->numIntVar; i++)
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
					temp->ubdsjncs++;
				}
				else if (config.BRN_STR == 1)
				{
					if (abs(fracVal - ceil(fracVal)) > 0.5)
					{
						temp->disjncs[i] = 1;
						temp->disjncsVal[i][0] = ceil(fracVal);
						temp->disjncsVal[i][1] = parent->disjncsVal[i][1];
						temp->lbdsjncs++;
					}
					else
					{
						temp->disjncs[i] = 1;
						temp->disjncsVal[i][0] = parent->disjncsVal[i][0];
						temp->disjncsVal[i][1] = floor(fracVal);
						temp->ubdsjncs++;
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
					temp->lbdsjncs++;
				}
				else if (config.BRN_STR == 1)
				{
					if (abs(fracVal - ceil(fracVal)) <= 0.5)
					{
						temp->disjncs[i] = 1;
						temp->disjncsVal[i][0] = ceil(fracVal);
						temp->disjncsVal[i][1] = parent->disjncsVal[i][1];
						temp->lbdsjncs++;
					}
					else
					{
						temp->disjncs[i] = 1;
						temp->disjncsVal[i][0] = parent->disjncsVal[i][0];
						temp->disjncsVal[i][1] = floor(fracVal);
						temp->ubdsjncs++;
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
	if (!(temp->vars = (dVector)arr_alloc(temp->numVar, double)))
		errMsg("allocation", "newNode", "temp->vars", 0);
	
	for (int v = 0; v < temp->numVar; v++)
		temp->vars[v] = rootNode->vars[v];

	return temp;

}//End newNode()

struct BnCnodeType *copyNode(struct BnCnodeType *node, double thresh)
{
	int i;
	struct BnCnodeType *temp = (struct BnCnodeType *)malloc(sizeof(struct BnCnodeType));

	temp->key = node->key;
	temp->parentkey = node->parentkey;
	temp->depth = node->depth;
	temp->isleft = node->isleft;
	temp->fracVal = node->fracVal;
	temp->ishrstic = true;
	temp->varId = node->varId;
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
int addBnCDisjnct(cellType *cell, dVector  *disjncsVal, struct BnCnodeType * node, dVector bdl, dVector bdu)
{
	int 	cnt;
	dVector	lbounds, ubounds;
	int numCols = node->numVar - 1; /* eta column excluded */

	if (!(lbounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "addDisjnct", "lbounds", 0);

	if (!(ubounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "addDisjnct", "ubounds", 0);

	/* Change the Bounds */
	for (cnt = 0; cnt < numCols; cnt++) {
		if (cnt >= node->stInt && cnt <= node->edInt) {
			lbounds[cnt] = max(bdl[cnt], disjncsVal[cnt][0] - config.INT_TOLERANCE);
			cell->master->bdl[cnt] = lbounds[cnt];
			ubounds[cnt] = min(bdu[cnt], disjncsVal[cnt][1] + config.INT_TOLERANCE);
			cell->master->bdu[cnt] = ubounds[cnt];
		}
		else {
			lbounds[cnt] = bdl[cnt];
			cell->master->bdl[cnt] = lbounds[cnt];
			ubounds[cnt] = bdu[cnt];
			cell->master->bdu[cnt] = ubounds[cnt];
		}

	}

	/* Adjust the initial incumbent of the node based on the disjunctions */
	if (!node->ishrstic && node->depth != 0)
	{
		if (node->isleft)
		{
			node->vars[node->varId + 1] = disjncsVal[node->varId][1];
			for (cnt = 0; cnt < numCols; cnt++) {
				if (node->disjncs[cnt] == 0) {
					node->vars[cnt+1] = 1;
					break;
				}
			}
		}
		else {
			node->vars[node->varId + 1] = disjncsVal[node->varId][0];
		}
	}

#if defined(BNC_CHECK)
	printf("\nlower bounds");
	for (int n = 0; n < node->numIntVar; n++)
		printf("%0.4f, ", cell->master->bdl[n]);
	printf("\nupper bounds");
	for (int n = 0; n < node->numIntVar; n++)
		printf("%0.4f, ", cell->master->bdu[n]);
#endif // defined(BNC_CHECK)


	if (changeQPbds(cell->master->lp, numCols, lbounds, ubounds, node->vars, 0)) {
		errMsg("algorithm", "algoIntSD", "failed to change the bounds to convert the problem into QP", 0);
		return 1;
	}

	mem_free(lbounds);
	mem_free(ubounds);

	return 0;
}

// Truncate the variables from SD based on the lower bounds and upper bounds of the original problem
void truncate(dVector var, dVector lb, dVector ub, int cnt)
{
	for (int v = 0; v < cnt; v++)
	{
		if (var[v+1] - config.INT_TOLERANCE <= lb[v]) var[v+1] = lb[v];
		if (var[v+1] + config.INT_TOLERANCE >= ub[v]) var[v+1] = ub[v];
	}
}

/* Return the previous active node of the input node 
   The reason for this subroutine is that some of the nodes are deactivated */
struct BnCnodeType *nextNode(struct BnCnodeType *node)
{
	struct BnCnodeType *temp = NULL;
	temp = node;

	while (true)
	{
		temp = temp->prevnode;

		if (temp->key == 0)
		{
			break;
		}

		if (temp->isActive == true)
		{
			break;
		}

	}


	return temp;
}

/* This function will reform the bnc cuts based on the most recent observations of omega passed in as _observ_, and the istar's
* which have already been stored in the _istar_ field of each cut. If an istar field does not exist for a given observation,
* it would be calculated using computeIstar(). */
void revisitNode(numType *num, coordType *coord, basisType *basis, sigmaType *sigma, deltaType *delta, omegaType *omega, sampleType *sample,
	dVector Xvect, int numSamples, bool *dualStableFlag, dVector pi_ratio, int numIter, double lb, oneCut *cut) {

	dVector 	piCbarX;
	double  cummOld = 0.0, cummAll = 0.0, argmax, multiplier;
	int	 	istar, idx, c, obs, sigmaIdx, lambdaIdx;
	bool    pi_eval_flag = false;


	/* Pre-compute pi x Cbar x x as it is independent of observations */
	if (!(piCbarX = arr_alloc(sigma->cnt, double)))
		errMsg("Allocation", "SDCut", "pi_Tbar_x", 0);
	for (c = 0; c < sigma->cnt; c++)
		piCbarX[c] = vXv(sigma->vals[c].piC, Xvect, coord->CCols, num->cntCcols);

	/* Adjust the alpha and beta of the input cut according to the most recent number of observations */
	cut->alpha *= cut->numSamples;

	for (c = 1; c <= num->prevCols; c++)
		cut->beta[c] *= cut->numSamples;

	/* Test for omega issues */
	for (obs = cut->omegaCnt; obs < omega->cnt; obs++) {
		/* identify the maximal Pi/basis that generates the maximal Pi for each observation */
		istar = computeIstar(num, coord, basis, sigma, delta, sample,
			piCbarX, Xvect, omega->vals[obs], obs, numSamples, pi_eval_flag, &argmax, false);

		if (istar < 0) {
			errMsg("algorithm", "SDCut", "failed to identify maximal Pi for an observation", 0);
			return;
		}
		cut->iStar[obs] = istar;

		if (num->rvdOmCnt > 0) {
			for (idx = 0; idx <= basis->vals[istar]->phiLength; idx++) {
				sigmaIdx = basis->vals[istar]->sigmaIdx[idx];
				lambdaIdx = sigma->lambdaIdx[sigmaIdx];
				if (idx == 0)
					multiplier = 1.0;
				else
					multiplier = omega->vals[obs][coord->rvOffset[2] + basis->vals[istar]->omegaIdx[idx]];

				/* Start with (Pi x bBar) + (Pi x bomega) + (Pi x Cbar) x X */
				cut->alpha += omega->weights[obs] * multiplier * (sigma->vals[sigmaIdx].pib + delta->vals[lambdaIdx][obs].pib);

				for (c = 1; c <= num->cntCcols; c++)
					cut->beta[coord->CCols[c]] += omega->weights[obs] * multiplier * sigma->vals[sigmaIdx].piC[c];
				for (c = 1; c <= num->rvCOmCnt; c++)
					cut->beta[coord->rvCOmCols[c]] += omega->weights[obs] * multiplier * delta->vals[lambdaIdx][obs].piC[c];
			}
		}
		else {
			cut->alpha += sigma->vals[istar].pib * omega->weights[obs];
			cut->alpha += delta->vals[sigma->lambdaIdx[istar]][obs].pib * omega->weights[obs];

			for (c = 1; c <= num->cntCcols; c++)
				cut->beta[coord->CCols[c]] += sigma->vals[istar].piC[c] * omega->weights[obs];
			for (c = 1; c <= num->rvCOmCnt; c++)
				cut->beta[coord->rvCols[c]] += delta->vals[sigma->lambdaIdx[istar]][obs].piC[c] * omega->weights[obs];
		}
	}


	cut->alpha = cut->alpha / numSamples;

	for (c = 1; c <= num->prevCols; c++)
		cut->beta[c] = cut->beta[c] / numSamples;
	cut->beta[0] = 1.0;			/* coefficient of eta coloumn */

	if(piCbarX) mem_free(piCbarX);


}//END reform_cuts

/* Return the previous node of the input node */
struct BnCnodeType *prevNode(int key)
{

	struct BnCnodeType *temp = rootNode;
	struct BnCnodeType *prevtemp = rootNode;

	while (true)
	{

		if (temp->nextnode == NULL || key == 0)
		{
			return NULL;
		}
		temp = temp->nextnode;

		if (temp->key == key)
		{
			return prevtemp;
		}
		prevtemp = prevtemp->nextnode;

	}

	return rootNode;
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
void fracLamda(cellType *cell, struct BnCnodeType *node) {
	int i, j, k, totLambda;

	node->Lambdasize = cell->lambda->cnt;
	totLambda = 0;

	for (j = 0; j < cell->activeCuts->vals[cell->iCutIdx]->numSamples; j++) {
		if (config.Pi_EVAL_FLAG == 1)
		{
			node->IncumbiStar[totLambda] = cell->activeCuts->vals[cell->iCutIdx]->iStar[j];
		}
		totLambda += 1;
	}


	/*
	* Find duplicate elements in array
	*/
	if (config.Pi_EVAL_FLAG == 1) {
		for (i = 0; i < totLambda - 1; i++)
		{
			for (j = i + 1; j < totLambda; j++)
			{
				/* If any duplicate found */
				if (node->IncumbiStar[i] == node->IncumbiStar[j])
				{
					/* Delete the current duplicate element */
					for (k = j; k < totLambda - 1; k++)
					{
						node->IncumbiStar[k] = node->IncumbiStar[k + 1];
					}

					/* Decrement size after removing duplicate element */
					totLambda--;

					/* If shifting of elements occur then don't increment j */
					j--;
				}
			}
		}
	}


	node->tightPi = totLambda < node->Lambdasize ? totLambda : node->Lambdasize;

	node->fracPi = ((double)totLambda) / ((double)node->Lambdasize);
}//END fracLamda()

/* sum the number of disjunctions if it is equal to the number of integer variables return one */
int sumDisjncs(iVector disjncs, int numIntVar) {
	for (int v = 0; v < numIntVar; v++) {
		if (disjncs[v] == 0) {
			return 1;
		}
	}

	return 0;
}// END sumDisjncs

/* Print the node information summary */
void printNodesummary(struct BnCnodeType *node) {
	printf("\n\n-----------------------------------------\n");
	printf("key: %d\n", node->key);
	printf("depth: %d\n", node->depth);
	printf("numVars: %d\n", node->numVar);
	printf("numIntvars: %d\n", node->numIntVar);
	printf("st int idx: %d\n", node->stInt);
	printf("ed int idx: %d\n", node->edInt);
	printf("disjunctions \n");
	printIntvec(node->disjncs, node->numIntVar, NULL);
	printf("imposed bounds: \n");
	for (int v = 0; v < node->numIntVar; v++) {
		printf("\nvar %d: ",v);
		printVector(node->disjncsVal[v], 2, NULL);
	}
	printf("vars: \n");
	printVector(node->vars, node->numVar, NULL);
	printf("\n\n-----------------------------------------\n");
}

/* Free the queue by passing the root node and using next nodes */
void freeNodes(struct BnCnodeType *root) {

	struct BnCnodeType *node = root;
	while (node != NULL) {
		struct BnCnodeType *next = node->nextnode;
		freeNode(node);
		node = next;
	}
	return;
}//End freeNode()

/* free a single node */
void freeNode(struct BnCnodeType *node) {

	// Return 1 when the tree is empty
	if (node == NULL) return;
	if (node->disjncs) mem_free(node->disjncs);
	if (node->disjncsVal) {
		for (int i = 0; i < node->numIntVar; i++)
			if (node->disjncsVal[i])  mem_free(node->disjncsVal[i]);
		mem_free(node->disjncsVal);
	}
	if (node->IncumbiStar) mem_free(node->IncumbiStar);
	if (node->ParIncumbiStar) mem_free(node->ParIncumbiStar);
	if (node->vars)  mem_free(node->vars);
	mem_free(node);

}//End freeNode()

 /* free a single node partially */
void freepartialNode(struct BnCnodeType *node) {

	node->isActive = false;

	// Return 1 when the tree is empty
	if (node == NULL) return;
	if (node->key == 0) return;
	if (node->disjncs) {
		mem_free(node->disjncs);
		node->disjncs = NULL;
	}
	if (node->disjncsVal) {
		for (int i = 0; i < node->numIntVar; i++)
			if (node->disjncsVal[i]) {
				mem_free(node->disjncsVal[i]);
				node->disjncsVal[i] = NULL;
			}
		mem_free(node->disjncsVal);
		node->disjncsVal = NULL;
	}
	if (node->IncumbiStar) {
		mem_free(node->IncumbiStar);
		node->IncumbiStar = NULL;
	}
	if (node->ParIncumbiStar) {
		mem_free(node->ParIncumbiStar);
		node->ParIncumbiStar = NULL;
	}


}//End freepartialNode()
