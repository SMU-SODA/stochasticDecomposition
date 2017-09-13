/*
 * twoSD.h
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#ifndef TWOSD_H_
#define TWOSD_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"
#include "stoc.h"
#include "cell.h"

#undef ALGO_CHECK

/* twoSD.c */
int solveSDCell(stocType *stoc, probType **prob, cellType *cell);
void writeSDStatistic(FILE *soln, probType **prob, cellType *cell, string probName, int numStages);

void verifySDSetup();

int solveSDMaster(numType *num, sparseVector *dBar, cellType *cell, int IniRow, double lb);
int changeEtaCol(LPptr lp, int numRows, int numCols, int k, cutsType *cuts, double lb);
int updateRHS(LPptr lp, cutsType *cuts, int numIter, double lb);
int changeEtaCol(LPptr lp, int numRows, int numCols, int k, cutsType *cuts, double lb);
int updateRHS(LPptr lp, cutsType *cuts, int numIter, double lb);
int formSDCut(probType *prob, cellType *cell, vector Xvect, int omegaIdx, BOOL newOmegaFlag, BOOL isIncumb);
oneCut *SDCut(numType *num, coordType *coord, basisType *basis, sigmaType *sigma, deltaType *delta, omegaType *omega, vector Xvect, int numSamples,
		BOOL *dualStableFlag, vector pi_ratio, double lb);
int checkImprovementSD(probType *prob, cellType *cell, int candidCut);

/* optimal.c */
BOOL optimalSD(probType **prob, cellType *cell);
BOOL preTest(cellType *cell);
BOOL fullTest(probType **prob, cellType *cell);
cutsType *chooseCuts(cutsType *cuts, vector pi, int lenX);
void reformCuts(basisType *basis, sigmaType *sigma, deltaType *delta, omegaType *omega, numType *num, coordType *coord, cutsType *gCuts, int *observ, int k, int lbType, int lb, int lenX);
double calcBootstrpLB(probType *prob, vector incumbX, vector piM, vector djM, int currIter, double quadScalar, cutsType *cuts);
void empiricalDistribution(omegaType *omega, int *cdf);
void resampleOmega(intvec cdf, intvec observ, int numSamples);

/* evaluate.c */
int evaluateSD(FILE **soln, stocType *stoc, probType **prob, cellType *cell, vector Xvect);

#endif /* TWOSD_H_ */
