/*
 * solver.h
 *
 *  Created on: Apr 20, 2014
 *      Author: gjharsha
 */
#ifndef MTSD_SOLVER_H_
#define MTSD_SOLVER_H_

#include <./ilcplex/cpxconst.h>
#include <utils.h>

#define		ENVptr			CPXENVptr
#define 	LPptr			CPXLPptr

#define		ON				CPX_ON
#define		OFF				CPX_OFF
#define 	INFBOUND    	CPX_INFBOUND

#define		PARAM_SCRIND	CPX_PARAM_SCRIND
#define		PARAM_SCAIND	CPX_PARAM_SCAIND
#define		PARAM_LPMETHOD	CPX_PARAM_LPMETHOD
#define		PARAM_QPMETHOD	CPX_PARAM_QPMETHOD
#define 	PARAM_PREIND	CPX_PARAM_PREIND

#define		ALG_AUTOMATIC	CPX_ALG_AUTOMATIC
#define		ALG_PRIMAL		CPX_ALG_PRIMAL
#define		ALG_DUAL		CPX_ALG_DUAL
#define		ALG_NET			CPX_ALG_NET
#define		ALG_BARRIER		CPX_ALG_BARRIER
#define		ALG_SIFTING		CPX_ALG_SIFTING
#define		ALG_CONCURRENT	CPX_ALG_CONCURRENT

#define		STAT_OPTIMAL	CPX_STAT_OPTIMAL
#define		STAT_INFEASIBLE	CPX_STAT_INFEASIBLE

#define		MSGBUFSIZE		CPXMESSAGEBUFSIZE

#define		PROB_LP			CPXPROB_LP
#define		PROB_QP			CPXPROB_QP
#define		PROB_MILP		CPXPROB_MILP
#define		PROB_MIQP		CPXPROB_MIQP

#define		AT_LOWER        CPX_AT_LOWER
#define		BASIC           CPX_BASIC
#define		AT_UPPER        CPX_AT_UPPER
#define		FREE_SUPER      CPX_FREE_SUPER

#define		MIP_OPTIMAL		CPXMIP_OPTIMAL
#define		MIP_OPTIMAL_TOL	CPXMIP_OPTIMAL_TOL
#define		MIP_INFEASIBLE	CPXMIP_INFEASIBLE
#define     MIP_OPTIMAL_TOL CPXMIP_OPTIMAL_TOL

#define 	THREADS			CPXPARAM_Threads

int solveProblem(LPptr lp, cString pname, int type, int *status);
int getProbType(LPptr lp);
double getObjective(LPptr lp, int type);
int getPrimal(LPptr lp, dVector X, int length);
double getPrimalPoint(LPptr lp, int idx);
int getDual(LPptr lp, dVector Pi, int length);
int getDualSlacks(LPptr lp, dVector Dj, int length);
int getBasis(LPptr lp, iVector cstat, iVector rstat);
int getBinvC(LPptr lp, int col, dVector a );
int changeCoef(LPptr lp, int row, int col, double val);
int changeObjx(LPptr lp, int cnt, iVector indices, dVector values);
int changeRHS(LPptr lp, int cnt, iVector indices, dVector values);
int changeBDS(LPptr lp, int cnt, iVector indices, cString lu, dVector bd);
int changeCol(LPptr lp, int column, dVector coef, int start, int stop);
int changeCtype(LPptr lp, int cnt, iVector indices, cString ctype);
int changeProbType(LPptr lp, int type);
int addRow(LPptr lp, int nzcnt, double inputRHS, char inputSense, int matbeg, iVector rmatind, dVector rmatval, cString rowname);
int addCol(LPptr lp, int nzcnt, double objx, int cmatbeg, iVector cmatind, dVector cmatval, double bdu, double bdl, cString colname);
int removeRow(LPptr lp, int begin, int end);

int createProblem(char *probname, LPptr *lp);
int readProblem(char *probpath, LPptr lp);
LPptr setupProblem(cString name, int type, int numcols, int numrows, int objsense, dVector objx, dVector rhsx, cString sense, iVector matbeg, iVector matcnt,
		iVector matind, dVector matval, dVector lb, dVector ub, dVector rngval, cString *colname, cString *rowname, cString ctype);
int loadProblem(CPXLPptr lp, int numcols, int numrows, int objsense, dVector objx, dVector rhsx, cString sense, iVector matbeg, iVector matcnt,
		iVector matind, dVector matval, dVector lb, dVector ub, dVector rngval);
int loadProbwNames(LPptr lp, int numcols, int numrows, int objsense, dVector objx, dVector rhsx, cString sense, iVector matbeg, iVector matcnt,
		iVector matind, dVector matval, dVector lb, dVector ub, dVector rngval, cString *colname, cString *rowname);
int copyQPseparable(LPptr lp, double *qsepvec);
LPptr cloneProblem(LPptr origLp);
int writeProblem(LPptr lp, char *filename);

void openSolver();
void closeSolver();
int setIntParam(int paramname, int paramvalue);
void solverErrmsg(int status);
int changeLPSolverType(int method);
int changeQPSolverType(int method);

int getProbName(LPptr lp, cString probName, int len);
int getObjSen(LPptr lp);
int getNumRows(LPptr lp);
int getNumCols(LPptr lp);
int getNumBinary(LPptr lp);
int getCtype(LPptr lp, int start, int end, cString ctype);
int getNumInt(LPptr lp);
int getNumnz(LPptr lp);
int getObjx(LPptr lp, int start, int end, dVector obj);
int getRhsx(LPptr lp, int start, int end, dVector rhs);
int getSense(LPptr lp, int start, int end, cString sense);
int getCols(LPptr lp, int start, int end, int *cmatbeg, iVector cmatind, dVector cmatval, int cmatspace);
int getLb(LPptr lp, int start, int end, dVector lb);
int getUb(LPptr lp, int start, int end, dVector ub);
int getObjName(LPptr lp, cString objname);
int getCstoreSize(LPptr lp, int start, int end);
int getColName(LPptr lp, int start, int end, cString *colname, cString colnamestore, int csize);
int getRstoreSize(LPptr lp, int start, int end);
int getRowName(LPptr lp, int start, int end, cString *rowname, cString rownamestore, int rsize);
int getBasisHead(LPptr lp, iVector head, dVector basicX);
int getBasisInvRow(LPptr lp, int i, dVector phi);
int getBasisInvCol(LPptr lp, int i, dVector phi);
int getBasisInvARow(LPptr lp, int i, dVector phi);
int getBasisInvACol(LPptr lp, int i, dVector phi);
int freeProblem(LPptr lp);

#endif /* MTSD_SOLVER_H_ */
