/*
 * cuts.c
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

int addCut2Master(cellType *cell, oneCut *cut, BOOL scaleCut, int lenX, double lb) {
	intvec 	indices;
	int 	cnt;

	if (!(indices = arr_alloc(lenX + 1, int)))
		errMsg("Allocation", "addcut2Master", "fail to allocate memory to coefficients of beta",0);
	for (cnt = 1; cnt <= lenX; cnt++)
		indices[cnt] = cnt - 1;
	indices[0] = lenX;

	if ( config.MASTER_TYPE == PROB_QP )
		cut->alphaIncumb = cut->alpha - vXv(cut->beta, cell->incumbX, NULL, lenX);

	/* check to see if there is room for the candidate cut, else drop a cut */
	if (cell->cuts->cnt == cell->maxCuts) {
		/* make room for the latest cut */
		if( reduceCuts(cell->master, cell->cuts, scaleCut, cell->candidX, cell->piM, lenX, lb, cell->k, &cell->iCutIdx, config.TOLERANCE) < 0 ) {
			errMsg("algorithm", "addCut2Master", "failed to add reduce cuts to make room for candidate cut", 0);
			return -1;
		}
	}

	/* add the cut to the cell cuts structure as well as on the solver */
	cell->cuts->vals[cell->cuts->cnt] = cut;
	if ( addRow(cell->master->lp, lenX + 1, cut->alphaIncumb, GE, 0, indices, cut->beta) ) {
		errMsg("solver", "addcut2Master", "failed to add new row to problem in solver", 0);
		return -1;
	}
	cut->rowNum = cell->master->mar++;

	mem_free(indices);
	return cell->cuts->cnt++;
}//END addCuts2Master()

/* This function loops through a set of cuts and find the highest cut height at the specified position x */
double maxCutHeight(cutsType *cuts, vector xk, int betaLen, BOOL scaleCut, int currIter, double lb) {
	double Sm = -INF, ht = 0.0;
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		ht = cutHeight(cuts->vals[cnt], xk, betaLen, scaleCut, currIter, lb);
		if (Sm < ht) {
			Sm = ht;
		}
	}

	return Sm;
}//END maxCutHeight

/* This function calculates and returns the height of a given cut at a given X.  It includes the k/(k-1) update, but does not include
 * the coefficients due to the cell. */
double cutHeight(oneCut *cut, vector xk, int betaLen, BOOL scaleCut, int currIter, double lb) {
	double height;

	/* A cut is calculated as alpha - beta x X */
	height = cut->alpha - vXv(cut->beta, xk, NULL, betaLen);

	if ( scaleCut ) {
		double t_over_k = ((double) cut->numSamples / (double) currIter);

		/* Weight cut based on number of observations used to form it */
		height *= t_over_k;

		/* Updated for optimality cut height*/
		height += (1 - t_over_k) * lb;
	}
	return height;
}//END cutHeight()

/* This function allocates memory for the arrays inside a single cut, and initializes its values accordingly.  The cut structure
 * itself is assumed to be already allocated.  Note, each beta vector contains room for its one-norm, thought it just gets filled
 * with zero anyway. */
oneCut *newCut(int numX, int numIstar, int numSamples) {
	oneCut *cut;

	cut = (oneCut *) mem_malloc (sizeof(oneCut));
	cut->numSamples = numSamples;
	cut->omegaCnt = numIstar;
	cut->isIncumb = FALSE; 								/* new cut is by default not an incumbent */
	cut->alphaIncumb = 0.0;
	cut->rowNum = -1;

	if (!(cut->iStar = arr_alloc(numIstar, int)))		/* when used in aggregate cut mode (MULTI_CUT = 0), this holds the index of agent cuts */
		errMsg("allocation", "new_cut", "iStar", 0);
	if (!(cut->beta = arr_alloc(numX + 1, double)))
		errMsg("allocation", "new_cut", "beta", 0);

	cut->alpha = 0.0;

	return cut;
}//END newCut

/* This function allocates memory for a new cut structure.  This entails the structure itself, and the _val_ array of oneCut pointers
 * inside the structure.  The actual oneCut structures are allocated according to the numBeta parameter, via calls to new_cut(). */
cutsType *newCuts(int maxCuts) {
	cutsType *cuts;

	if (maxCuts == 0)
		return NULL;

	if (!(cuts = (cutsType *) mem_malloc (sizeof(cutsType))))
		errMsg("allocation", "newCuts", "cuts",0);
	if (!(cuts->vals = (oneCut **) arr_alloc (maxCuts, oneCut)))
		errMsg("allocation", "newCuts", "oneCuts",0);
	cuts->cnt = 0;

	return cuts;
}//END newCuts

/* This function will remove the oldest cut whose corresponding dual variable is zero (thus, a cut which was slack in last solution). */
int reduceCuts(oneProblem *master, cutsType *cuts, BOOL scaleCut, vector candidX, vector pi, int betaLen, double lb, int currentIter, int *iCutIdx, double TOLERANCE) {
	double height, minHeight;
	int minObs, oldestCut,idx;

	minObs 	  = currentIter;
	oldestCut = cuts->cnt;

	/* identify the oldest loose cut */
	for (idx = 0; idx < cuts->cnt; idx++) {
		if ( idx == (*iCutIdx) || cuts->vals[idx]->rowNum < 0)
			/* avoid dropping incumbent cut and newly added cuts */
			continue;

		if (cuts->vals[idx]->numSamples < minObs && DBL_ABS(pi[cuts->vals[idx]->rowNum + 1]) <= TOLERANCE ) {
			minObs = cuts->vals[idx]->numSamples;
			oldestCut = idx;
		}
	}

	/* if the oldest loose cut is the most recently added cut, then the cut with minimium cut height will be dropped */
	if ( oldestCut == cuts->cnt ) {
		//minHeight = cutHeight(lbType, cell[agentIdx]->cuts->vals[0], cell[agentIdx]->k, candidX, betaLen, lb);
		minHeight = cutHeight(cuts->vals[0], candidX, betaLen, scaleCut, currentIter, lb);
		oldestCut = 0;

		for (idx = 1; idx < cuts->cnt; idx++) {
			if (idx == (*iCutIdx))
				continue;

			//height = cutHeight(lbType, cell[agentIdx]->cuts->vals[idx], cell[agentIdx]->k, candidX, betaLen, lb);
			height = cutHeight(cuts->vals[idx], candidX, betaLen, scaleCut, currentIter, lb);
			if (height < minHeight) {
				minHeight = height;
				oldestCut = idx;
			}
		}
	}

	/* drop the selected cut and swap the last cut into its place */
	if ( dropCut(master, cuts, oldestCut, iCutIdx) ){
		errMsg("algorithm", "reduceCuts", "failed to drop a cut", 0);
		return -1;
	}

	return oldestCut;
}//END reduceCuts()

/* This function removes a cut from both the cutType structure and the master problem constraint matrix.  In the cuts->vals array, the last
 * cut is swapped into the place of the exiting cut.  In the constraint matrix, the row is deleted, and the row numbers of all constraints
 * below it are decremented. */
int dropCut(oneProblem *master, cutsType *cuts, int cutIdx, int *iCutIdx) {
	int idx, deletedRow;

	deletedRow = cuts->vals[cutIdx]->rowNum;
	/* Get rid of the indexed cut on the solver */
	if (  removeRow(master->lp, deletedRow, deletedRow) ) {
		errMsg("solver", "dropCut", "failed to remove a row from master problem", 0);
		return 1;
	}
	freeOneCut(cuts->vals[cutIdx]);

	/* move the last cut to the deleted cut's position (structure) */
	cuts->vals[cutIdx] = cuts->vals[--cuts->cnt];

	/* if the swapped cut happens to be the incumbent cut, then update its index */
	if ( (*iCutIdx) == cuts->cnt )
		(*iCutIdx) = cutIdx;

	for (idx = 0; idx < cuts->cnt; idx++) {
		if (cuts->vals[idx]->rowNum > deletedRow)
			--cuts->vals[idx]->rowNum;
	}

	/* decrease the number of rows on solver */
	master->mar--;

	return 0;
}//END dropCut()

void freeOneCut(oneCut *cut) {

	if (cut) {
		if (cut->iStar)
			mem_free(cut->iStar);
		if (cut->beta)
			mem_free(cut->beta);
		mem_free(cut);
	}
}

void freeCutsType(cutsType *cuts) {
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++)
		freeOneCut(cuts->vals[cnt]);
	mem_free(cuts->vals);
	mem_free(cuts);
}//END freeCuts
