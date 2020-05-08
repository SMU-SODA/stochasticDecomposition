/*
 * solver.c
 *
 *  Created on: Apr 20, 2014
 *      Author: Harsha Gangammanavar
 */

#include <utils.h>
#include <solver.h>

extern cString 	outputDir;
ENVptr	env;

int solveProblem(LPptr lp, cString pname, int type, int *status, double mipGap) {
	int		aggres = 0, res = 0;

	solveagain:
	switch  ( type ) {
	case PROB_LP:
		setIntParam(PARAM_PREIND, OFF);
		(*status) = CPXlpopt(env, lp);
		setIntParam(PARAM_PREIND, ON);
		break;
	case PROB_QP:
		(*status) = CPXqpopt(env, lp);
		break;
	case PROB_MILP:
		setDoubleParam(CPXPARAM_MIP_Tolerances_MIPGap, mipGap);
		(*status) = CPXmipopt(env, lp);
		break;
	case PROB_MIQP:
		(*status) = CPXmipopt(env, lp);
		break;
	default:
		break;
	}

	RESOLVE:
	(*status) = CPXgetstat(env, lp);
	if ((*status) != STAT_OPTIMAL && (*status) != MIP_OPTIMAL && (*status)!= MIP_OPTIMAL_TOL  ) {
		writeProblem(lp, "error.lp");
		if ((*status) == STAT_INFEASIBLE || (*status) == MIP_INFEASIBLE) {
			fprintf(stderr, "\nWarning: Problem is infeasible. ");
			return (*status);
		}
		else if ( (type == PROB_LP || type == PROB_QP) && (*status) == 6 ){
			setIntParam(PARAM_SCAIND, 1);
			aggres++;
			if(aggres == 1)
				goto solveagain;
			else
				goto skip;
		} else if ( (type == PROB_QP) && (*status) == 2 ) {
			changeLPSolverType(ALG_PRIMAL);
			(*status) = CPXlpopt(env, lp);
			goto RESOLVE;
		} else if ( type == PROB_LP && (*status) == 12 ) {
			switch (res) {
			case 0: changeLPSolverType(ALG_PRIMAL); res++; break;
			case 1: changeLPSolverType(ALG_PRIMAL); res++; break;
			case 2: changeLPSolverType(ALG_BARRIER); res++; break;
			default: fprintf(stderr, "Failed to resolve conflict.\n"); return (*status);
			}
			(*status) = CPXlpopt(env, lp);
			goto RESOLVE;
		} else if ( type == PROB_LP && (*status) == 5 ) {
			setDoubleParam(CPX_PARAM_EPRHS, 1e-09); 		// Decrease the feasibility tolerance
			(*status) = CPXlpopt(env, lp);
			setDoubleParam(CPX_PARAM_EPRHS, 1e-06);	 		// Set it back to the default value 1e-06
			goto RESOLVE;
		}
		else {
			solverErrmsg((*status));
			return (*status);
		}
	}

	skip: if (aggres)
		setIntParam(PARAM_SCAIND, 0);

	return 0;
}//END solveProblem()

int solveProblemCallback(LPptr lp, cString pname, int *status, double mipGap) {
	int		aggres = 0, res = 0;

	setDoubleParam(CPXPARAM_MIP_Tolerances_MIPGap, mipGap);
	(*status) = CPXmipopt(env, lp);

	(*status) = CPXgetstat(env, lp);
	if ( (*status) != MIP_OPTIMAL && (*status)!= MIP_OPTIMAL_TOL ) {
		writeProblem(lp, "error.lp");
		if ((*status) == STAT_INFEASIBLE || (*status) == MIP_INFEASIBLE) {
			fprintf(stderr, "\nProblem is infeasible.\n");
			return (*status);
		}
		else {
			solverErrmsg((*status));
			return (*status);
		}
	}

	return 0;
}//END solveProblem()

int getProbType(LPptr lp) {
	int status;

	status = CPXgetprobtype (env, lp);
	if ( status < -1 ) {
		solverErrmsg(status);
		exit(1);
	}
	else
		return status;
}//END

double getObjective(LPptr lp, int type) {
	double ans;
	int status;

	/* Get the objective function value only */
	CPXsolution(env, lp, &status, &ans, NULL, NULL, NULL, NULL);
	if ( type == PROB_LP || type == PROB_QP ) {
		if(status != STAT_OPTIMAL && status != 6) {
			solverErrmsg(status);
			exit(1);
		}
	}
	else
		if(status != MIP_OPTIMAL && status != MIP_OPTIMAL_TOL ) {
			solverErrmsg(status);
			exit(1);
		}

	return ans;
}//END getObjective()

int getPrimal(LPptr lp, dVector X, int length) {
	int status;

	status = CPXgetx(env, lp, X+1, 0, length-1);
	if ( status )
		solverErrmsg(status);
	else
		X[0] = oneNorm(X+1, length);

	return status;
}//END getPrimal()

double getPrimalPoint(LPptr lp, int idx) {
	int 	status;
	double 	x[1];

	status = CPXgetx(env, lp, x, idx, idx);
	if ( status )
		solverErrmsg(status);

	return x[0];
}//END getPrimalPoint

int getDual(LPptr lp, dVector Pi, int length) {
	int	status;

	status = CPXgetpi(env, lp, Pi+1, 0, length-1);
	if(status)
		solverErrmsg(status);
	else
		Pi[0] = oneNorm(Pi+1, length);

	return status;
}//END getDual()

int getDualSlacks(LPptr lp, dVector Dj, int length){
	int	status;

	status = CPXgetdj(env, lp, Dj+1, 0, length-1);
	if(status)
		solverErrmsg(status);
	else
		Dj[0] = oneNorm(Dj+1, length);

	return status;
}//END getDualSlacks()

int getBasis(LPptr lp, iVector cstat, iVector rstat){
	int status;

	if ( rstat != NULL && cstat != NULL ) {
		status = CPXgetbase(env, lp, cstat+1, rstat+1);
		if ( status )
			solverErrmsg(status);
	}
	else if ( rstat != NULL ) {
		status = CPXgetbase(env, lp, NULL, rstat+1);
		if ( status )
			solverErrmsg(status);
	}
	else if ( cstat != NULL ) {
		status = CPXgetbase(env, lp, cstat+1, NULL);
		if ( status )
			solverErrmsg(status);
	}

	return status;
}//END get_basis()

void openSolver(){
	int 	status = 0;

	env = CPXopenCPLEX(&status);
	if ( env == NULL ) {
		solverErrmsg(status);
		errMsg("solver", "openSolver", "could not open CPLEX environment", 0);
		goto TERMINATE;
	}

	status = setIntParam(PARAM_SCRIND, OFF);
	if ( status ) {
		errMsg("solver", "open_solver", "screen output", 0);
		goto TERMINATE;
	}

	TERMINATE:
	if ( status )
		closeSolver();

}// openSolver()

void closeSolver(){
	int status;

	if (env != NULL) {
		status = CPXcloseCPLEX(&env);
		if ( status ) {
			solverErrmsg(status);
			errMsg("solver", "close_solver", "could not close CPLEX environment", 1);
		}
	}

}//END closeSolver()

int changeLPSolverType(int method) {
	int status = 0;

	status = setIntParam(PARAM_LPMETHOD, method);
	if (status) {
		solverErrmsg(status);
		return 1;
	}
	return 0;
}//END changeLPSolverType()

int changeQPSolverType(int method) {
	int status = 0;

	status = setIntParam(PARAM_QPMETHOD, method);
	if (status) {
		solverErrmsg(status);
		return 1;
	}
	return 0;
}//END changeQPSolverType()

int setDoubleParam (int paramName, double paramVal) {
	return CPXsetdblparam(env, paramName, paramVal);
}

int setIntParam(int paramName, int paramValue) {
	return CPXsetintparam (env, paramName, paramValue);
}//END setIntParam()

int getIntParam(int paramName, int *paramValue) {
	return CPXgetintparam (env, paramName, paramValue);
}//END setIntParam()

void solverErrmsg(int status){
	char  errmsg[CPXMESSAGEBUFSIZE];

	CPXgeterrorstring (env, status, errmsg);
	fprintf (stderr, "%s", errmsg);

}//END solverErrmsg()

int createProblem(char *probname, LPptr *lp) {
	int status;

	(*lp) = CPXcreateprob (env, &status, probname);
	if ( (*lp) == NULL )
		solverErrmsg(status);

	return status;
}// END createProb()

int readProblem(cString probpath, LPptr lp) {
	int status;

	status = CPXreadcopyprob (env, lp, probpath, "MPS");
	if( status ) {
		status = CPXreadcopyprob(env, lp, probpath,"LP");
		if ( status )
			solverErrmsg(status);
	}

	return status;
}// END readProb()

LPptr setupProblem(cString name, int type, int numcols, int numrows, int objsense, dVector objx, dVector rhsx, cString sense, iVector matbeg, iVector matcnt,
		iVector matind, dVector matval, dVector lb, dVector ub, dVector rngval, cString *colname, cString *rowname, cString ctype) {
	LPptr	lp = NULL;
	int		*indices, c;

	createProblem(name, &lp);
	if ( lp == NULL ) {
		errMsg("solver", "setupProblem", name, 1);
		return NULL;
	}

	if ( loadProbwNames(lp, numcols, numrows, objsense, objx, rhsx, sense, matbeg, matcnt, matind, matval, lb, ub, rngval, colname, rowname) ) {
		errMsg("solver", "setupProblem", "failed to load problem onto the solver", 0);
		return NULL;
	}


	if ( type == PROB_MILP || type == PROB_MIQP ) {
		if ( !(indices = (iVector) arr_alloc(numcols, int)) )
			errMsg("allocation", "setupProblem", "indices", 0);
		for ( c = 0; c < numcols; c++ )
			indices[c] = c;
		if ( changeCtype(lp, numcols, indices, ctype) ) {
			errMsg("solver", "setupProblem", "failed to change column type", 0);
			return NULL;
		}
		mem_free(indices);
	}


	if ( changeProbType(lp, type) ) {
		errMsg("Problem Setup", "new_subprob", "subprob",0);
		return NULL;
	}


	return lp;
}//END setupProblem()

int loadProblem(LPptr lp, int numcols, int numrows, int objsense, dVector objx, dVector rhsx, cString sense, iVector matbeg, iVector matcnt,
		iVector matind, dVector matval, dVector lb, dVector ub, dVector rngval) {
	int status;
	status = CPXcopylp(env, lp, numcols, numrows, objsense, objx, rhsx, sense, matbeg, matcnt, matind, matval, lb, ub, rngval);
	if(status)
		solverErrmsg(status);

	return status;
}// END loadwnames_prob()

int loadProbwNames(LPptr lp, int numcols, int numrows, int objsense, dVector objx, dVector rhsx, cString sense, iVector matbeg, iVector matcnt,
		iVector matind, dVector matval, dVector lb, dVector ub, dVector rngval, cString *colname, cString *rowname) {
	int status;

	status = CPXcopylpwnames (env, lp, numcols, numrows, objsense, objx, rhsx, sense, matbeg, matcnt, matind, matval, lb, ub, rngval, colname, rowname);
	if(status)
		solverErrmsg(status);

	return status;
}// END loadwnames_prob()

int writeProblem(LPptr lp, cString filename) {
	int status;
	char buffer[2*BLOCKSIZE];

	strcpy(buffer,outputDir);
	strcat(buffer,filename);

	status = CPXwriteprob (env, lp, buffer, "LP");
	if (status)
		solverErrmsg(status);

	return status;
}// END write_prob()

LPptr cloneProblem(LPptr origLp) {
	int status = 0;
	LPptr cloneLp = NULL;

	cloneLp = CPXcloneprob(env, origLp, &status);
	if (status != 0 || cloneLp == NULL)
		solverErrmsg(status);

	return cloneLp;
}//END cloneProb()

int getProbName(LPptr lp, cString probName, int len) {
	int status, surplus;

	status = CPXgetprobname (env, lp, probName, len, &surplus);
	if ( surplus < 0)
		solverErrmsg(status);

	return status;
}//END getProbName()

int getObjSen(LPptr lp) {
	int status;

	status = CPXgetobjsen(env, lp);
	if ( !(status) )
		solverErrmsg(status);

	return status;
}// END getObjsen()

int getNumRows(LPptr lp) {
	int status;

	status = CPXgetnumrows(env,lp);
	if ( !(status) )
		solverErrmsg(status);

	return status;
}// END getNumRows()

int getNumCols(LPptr lp) {
	int status;

	status = CPXgetnumcols(env, lp);
	if ( !(status) )
		solverErrmsg(status);

	return status;
}// END getNumCols()

int getNumBinary(LPptr lp){

	return CPXgetnumbin(env, lp);

}//END getNumBinary()

int getNumInt(LPptr lp){

	return CPXgetnumint(env, lp);

}//END getNumBinary()

int getCtype(LPptr lp, int start, int end, cString ctype) {
	int status;

	status = CPXgetctype (env, lp, ctype, start, end-1);
	if ( status )
		solverErrmsg(status);

	return status;
}//END getCtype()

int getNumnz(LPptr lp) {

	return CPXgetnumnz (env, lp);

}//END getNumnz

int getObjx(LPptr lp, int start, int end, dVector obj) {
	int status;

	if( (status = CPXgetobj (env, lp, obj, start, end-1)) )
		solverErrmsg(status);

	return status;
}// END getObjx()

int getRhsx(LPptr lp, int start, int end, dVector rhs) {
	int status;

	if( (status = CPXgetrhs (env, lp, rhs, start, end-1)))
		solverErrmsg(status);

	return status;
}// END getRhsx()

int getSense(LPptr lp, int start, int end, cString sense) {
	int status;

	if ( (status = CPXgetsense (env, lp, sense, start, end-1)))
		solverErrmsg(status);

	return status;
}// END getSense()

int getCols(LPptr lp, int start, int end, iVector cmatbeg, iVector cmatind, dVector cmatval, int cmatspace) {
	int status, nzcnt, surplus;

	status = CPXgetcols (env, lp, &nzcnt, cmatbeg, cmatind, cmatval, cmatspace, &surplus, start, end-1);
	if ( status != 0 || surplus < 0 ) {
		fprintf(stderr, "Surplus = %d", surplus);
		return 1;
	}
	if (nzcnt > cmatspace)
		printf("The allocated size and used size for constraint elements do not match.\n");

	return 0;
}// END getCols()

/* function reads the lower bound on variables indexed from _start_ to _end_ of problem input and stores them in _lb_. Returns 0 on success and 1 on
 * failure. */
int getLb(LPptr lp, int start, int end, dVector lb) {
	int status;

	if ( (status = CPXgetlb(env, lp, lb, start, end - 1)) )
		solverErrmsg(status);

	return status;
}//END getLb()

/* function reads the upper bound on variables indexed from _start_ to _end_ of problem input and stores them in _lb_. Returns 0 on success and 1 on
 * failure. */
int getUb(LPptr lp, int start, int end, dVector ub) {
	int status;

	if ( (status = CPXgetub(env, lp, ub, start, end - 1)))
		solverErrmsg(status);

	return status;
}//END getUb()

/* read problem objective function name into the cString _objname. Returns 0 on success and 1 on failure.*/
int getObjName(LPptr lp, cString objname) {

	return CPXgetobjname (env, lp, objname, NAMESIZE, NULL);

}// END getObjName()

int getCstoreSize(LPptr lp, int start, int end) {
	int status, surplus;

	status = CPXgetcolname (env, lp, NULL, NULL, 0, &surplus, 0, end - 1);
	if (( status != CPXERR_NEGATIVE_SURPLUS ) && ( status != 0 ))
		return 1;

	return surplus;
}// END getCstoreSize()

int getColName(LPptr lp, int start, int end, cString *colname, cString colnamestore, int csize) {
	int status;

	return CPXgetcolname (env, lp, colname, colnamestore, csize, &status, start, end-1);

}// END get_colname()

int getRstoreSize(LPptr lp, int start, int end) {
	int status, surplus;

	status = CPXgetrowname (env, lp, NULL, NULL, 0, &surplus, 0, end - 1);
	if (( status != CPXERR_NEGATIVE_SURPLUS ) && ( status != 0 ))
		return 1;

	return surplus;
}//END getRstoreSize()

int getRowName(LPptr lp, int start, int end, cString *rowname, cString rownamestore, int rsize) {
	int status;

	return CPXgetrowname (env, lp, rowname, rownamestore, rsize, &status, start, end-1);

}// END getRowName()

int getBasisHead(LPptr lp, iVector head, dVector basicX) {
	int status;

	status = CPXgetbhead(env, lp, head, basicX);
	if ( status )
		solverErrmsg(status);

	return status;
}//END getBasisHead()

int getBasisInvRow(LPptr lp, int i, dVector phi) {
	int status;

	status = CPXbinvrow(env, lp, i, phi);
	if ( status )
		solverErrmsg(status);

	return status;
}//END getBasicInvRow()

int getBasisInvCol(LPptr lp, int i, dVector phi) {
	int status;

	status = CPXbinvcol(env, lp, i, phi);
	if ( status )
		solverErrmsg(status);

	return status;
}//END getBasicInvRow()

int getBasisInvARow(LPptr lp, int i, dVector phi) {
	int status;

	status = CPXbinvarow(env, lp, i, phi);
	if ( status )
		solverErrmsg(status);

	return status;
}//END getBasicInvRow()

int getBasisInvACol(LPptr lp, int i, dVector phi) {
	int status;

	status = CPXbinvacol(env, lp, i, phi);
	if ( status )
		solverErrmsg(status);

	return status;
}//END getBasicInvRow()

int copyQPseparable(LPptr lp, dVector qsepvec){
	int status;
	/* NOTE: CPLEX evaluates the corresponding objective with a factor of 0.5 in front of the quadratic objective term.*/
	status = CPXcopyqpsep (env, lp, qsepvec);
	if(status)
		solverErrmsg(status);

	return status;
}//END copy_qp_separable()

int changeCoef(LPptr lp, int row, int col, double val) {
	int status;

	status = CPXchgcoef (env, lp, row, col, val);
	if ( status )
		solverErrmsg(status);

	return status;
}//END changeObjx()

int changeObjx(LPptr lp, int cnt, iVector indices, dVector values){
	int status;

	status = CPXchgobj (env, lp, cnt, indices, values);
	if ( status )
		solverErrmsg(status);

	return status;
}//END changeObjx()

int changeRHS(LPptr lp, int cnt, iVector indices, dVector values) {
	int status;

	status = CPXchgrhs (env, lp, cnt, indices, values);
	if(status)
		solverErrmsg(status);

	return status;
}//END changeRHS

int changeBDS(LPptr lp, int cnt, iVector indices, cString lu, dVector bd) {
	int status;

	status = CPXchgbds (env, lp, cnt, indices, lu, bd);
	if( status )
		solverErrmsg(status);

	return status;
}//END changeBDS

/* Change constraint sense solver function */
int changeSense(LPptr lp, int cnt, iVector indices, cString sense) {
	int status;

	status = CPXchgsense(env, lp, cnt, indices, sense);
	if (status)
		solverErrmsg(status);

	return status;
}//END changeSense()

int changeCol(LPptr lp, int column, dVector coef, int start, int stop){
	int		row, status;

	for (row = start; row < stop; row++){
		status = CPXchgcoef(env, lp, row, column, coef[row-start]);
		if( status )
			return status;
	}

	return 0;
}//END changeCol()

int changeCtype(LPptr lp, int cnt, iVector indices, cString ctype){
	int status;

	status = CPXchgctype (env, lp, cnt, indices, ctype);
	if ( status )
		solverErrmsg(status);

	return status;
}//END changeCtype()

int changeProbType(LPptr lp, int type) {
	int status;

	status = CPXchgprobtype (env, lp, type);
	if ( status )
		solverErrmsg(status);

	return status;
}//END changeProbType()

int addRow(LPptr lp, int nzcnt, double inputRHS, char inputSense, int matbeg, iVector rmatind, dVector rmatval,
		cString rowname) {
	cString *rNames;
	char	sense[1] = {'G'};
	double	rhs[1];
	int		status, rmatbeg[1];

	/* Row name */
	rNames = (cString *) arr_alloc(1, cString);
	rNames[0] = (cString) arr_alloc(NAMESIZE, char);
	strcpy(rNames[0],rowname);

	/* Row parameters */
	rhs[0] = inputRHS;
	sense[0] = inputSense;
	rmatbeg[0] = matbeg;

	/* Add row to the solver */
	status = CPXaddrows (env, lp, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rNames);
	if ( status )
		solverErrmsg(status);

	mem_free(rNames[0]); mem_free(rNames);

	return status;
}//END addRow()

int addCol(LPptr lp, int nzcnt, double objx, int matbeg, iVector cmatind, dVector cmatval, double bdu,
		double bdl, cString colname){
	cString	*cNames;
	int 	status, cmatbeg[1];
	double 	obj[1], lb[1], ub[1];

	/* Column name */
	cNames = (cString *) arr_alloc(1, cString);
	cNames[0] = (cString) arr_alloc(NAMESIZE, char);
	strcpy(cNames[0], colname);

	/* Column parameters */
	cmatbeg[0] = matbeg;
	obj[0] = objx;
	ub[0] = bdu;
	lb[0] = bdl;

	/* Add column to the solver */
	status = CPXaddcols(env, lp, 1, nzcnt, obj, cmatbeg, cmatind, cmatval, lb, ub, cNames);
	if ( status )
		solverErrmsg(status);

	mem_free(cNames[0]); mem_free(cNames);

	return status;
}//END addCol()

int removeRow(LPptr lp, int begin, int end){
	int status;

	status = CPXdelrows(env, lp, begin, end);
	if ( status )
		solverErrmsg(status);

	return status;
}//END removeRow

int freeProblem(LPptr lp) {
	int status;

	status = CPXfreeprob(env, &lp);
	if (status)
		solverErrmsg(status);

	return status;
}//END freeProb()

int refineConflict (LPptr lp , int mar, int mac) {
	int confnumrows_p, confnumcols_p, confstat_p;
	iVector rowind,rowbdstat, colind, colbdstat;

	if ( CPXrefineconflict(env, lp, &confnumrows_p, &confnumcols_p ) ) {
		fprintf (stderr, "Failed to obtain refine conflict info.\n");
		return 1;
	}

	/* Allocate memory to elements necessary for diagnosis. */
	rowind=(iVector)arr_alloc(mar, int);
	rowbdstat=(iVector)arr_alloc(mar, int);
	colind=(iVector)arr_alloc(mac, int);
	colbdstat=(iVector)arr_alloc(mac, int);

	if ( CPXXgetconflict(env, lp, &confstat_p, rowind, rowbdstat, &confnumrows_p, colind, colbdstat, &confnumcols_p) ) {
		fprintf (stderr, "Failed to obtain get conflict info.\n");
		return 1;
	}

	/* Identifying the type of conflict */
	switch (confstat_p){
	case CPX_STAT_CONFLICT_MINIMAL:
		printf("\tConflict is minimal\n");
		break;
	case CPX_STAT_CONFLICT_FEASIBLE:
		printf("\tNo conflict is available\n");
		break;
	case CPX_STAT_CONFLICT_ABORT_IT_LIM:
		printf("\tConflict resolution aborted\n");
		break;
	default :
		printf("\tUnknown Conflict Status\n");
		return 1;
	}

	/* Identifying type of row conflict */
	for (int n=0; n< confnumrows_p ;n++){
		switch (rowbdstat[n]){
		case CPX_CONFLICT_MEMBER:
			printf("\t\tRow id = %d,  Row status = conflict member\n",rowind[n]);
			break;
		case CPX_CONFLICT_LB:
			printf("\t\tRow id = %d,  Row status = conflict LB\n",rowind[n]);
			break;

		case CPX_CONFLICT_UB:
			printf("\t\tRow id = %d,  Row status = conflict UB\n",rowind[n]);
			break;
		case CPX_CONFLICT_POSSIBLE_MEMBER:
			printf("\t\tRow id = %d,  Row status = conflict possible member\n",rowind[n]);
			break;
		case CPX_CONFLICT_POSSIBLE_LB:
			printf("\t\tRow id = %d,  Row status = conflict possible LB\n",rowind[n]);
			break;
		case CPX_CONFLICT_POSSIBLE_UB:
			printf("\t\tRow id = %d,  Row status = conflict possible UB\n",rowind[n]);
			break;
		default :
			printf("\t\tUnknown row conflict status\n");
			return 1;
		}

	}

	/* Identifying the column conflicts */
	for (int n=0; n< confnumcols_p ;n++){
		switch (rowbdstat[n]){
		case CPX_CONFLICT_MEMBER:
			printf("\t\tColumn id = %d,  Column status = conflict member\n",colind[n]);
			break;
		case CPX_CONFLICT_LB:
			printf("\t\tColumn id = %d,  Column status = conflict LB\n",colind[n]);
			break;

		case CPX_CONFLICT_UB:
			printf("\t\tColumn id = %d,  Column status = conflict UB\n",colind[n]);
			break;
		case CPX_CONFLICT_POSSIBLE_MEMBER:
			printf("\t\tColumn id = %d,  Column status = conflict possible member\n",colind[n]);
			break;
		case CPX_CONFLICT_POSSIBLE_LB:
			printf("\t\tColumn id = %d,  Column status = conflict possible LB\n",colind[n]);
			break;

		case CPX_CONFLICT_POSSIBLE_UB:
			printf("\t\tColumn id = %d,  Column status = conflict possible UB\n",colind[n]);
			break;
		default :
			printf("\t\tUnknown columns conflict status\n");
			break;
		}
	}

	mem_free(rowind);
	mem_free(rowbdstat);
	mem_free(colind);
	mem_free(colbdstat);

	return 0;
}//END refineConflict()

/* Callback functions */

int setsolvecallbackfunc (void *solvecallback, void *cbhandle) {

	return CPXsetsolvecallbackfunc(env, solvecallback, cbhandle);

}//END setsolvecallbackfunc()

int getcallbacknodelp (void * cbdata, int wherefrom, CPXLPptr * nodelp) {

	return CPXgetcallbacknodelp (env, cbdata, wherefrom, nodelp);

}//END getcallbacknodelp()

int getcallbackinfo(void * cbdata, int wherefrom, int whichinfo, void * result_p) {

	return CPXgetcallbackinfo(env, cbdata, wherefrom, whichinfo, result_p);

}//END getcallbackinfo()

int getCallbackPrimal(void * cbdata, int wherefrom, dVector X, int length) {
	int status;
	status = CPXgetcallbacknodex (env, cbdata, wherefrom, X+1, 0, length-1);
	if ( status )
		solverErrmsg(status);
	else
		X[0] = oneNorm(X+1, length);

	return status;
}//END getCallbackPrimal()
