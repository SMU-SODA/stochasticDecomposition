/*
 * algo.h
 *
 *  Created on: Sep 29, 2015
 *      Author: Harsha Gangammanavar
 */

#ifndef BENDERS_H_
#define BENDERS_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"
#include "stoc.h"
#include "cell.h"

#undef ALGO_CHECK
#define DETAILED

/* subroutines in benders.c */
int solveBendersCell(stocType *stoc, probType **prob, cellType *cell);
void updateOmega(stocType *stoc, omegaType *omega);
void writeBendersStatistic(FILE *soln, probType **prob, cellType *cell);

int solveBendersMaster(numType *num, sparseVector *dBar, cellType *cell);
int checkImprovementBenders(probType *prob, cellType *cell, int candidCut);
int addCut2BendersMaster(cellType *cell, oneCut *cut, int lenX);
int formBendersCutPct(probType **prob, cellType *cell, vector Xvect, BOOL isIncumb);
int formBendersCutCnt(probType **prob, cellType *cell, vector Xvect, BOOL isIncumb);

int resolveInfeasibility(probType **prob, cellType *cell, BOOL newOmegaFlag, int omegaIdx);

/* optimal.c */
BOOL optimalBenders(probType **prob, cellType *cell);

/* evaluate.c */
int evaluateBenders(FILE **soln, stocType *stoc, probType **prob, cellType *cell, vector Xvect);

#endif /* BENDERS_H_ */
