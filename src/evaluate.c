/*
 * evaluate.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"

extern configType config;

int evaluate(FILE *soln, stocType *stoc, probType **prob, oneProblem *subprob, dVector Xvect) {
	dVector 	observ, rhs, cost, costTemp;
	iVector   objxIdx;
	double 	obj, mean, variance, stdev, temp;
	int		cnt, status, m;

	if ( !(observ = (dVector) arr_alloc(stoc->numOmega + 1, double)) )
		errMsg("allocation", "evaluate", "observ", 0);

	printf("\nStarting evaluation.\n");

	/* initialize parameters used for evaluations */
	cnt = 0.0; mean = 0.0; variance = 0.0; stdev = INFBOUND; cnt = 0;

	/* right-hand side */
	if (!(rhs =(dVector) arr_alloc(prob[1]->num->rows+1, double)))
		errMsg("Allocation", "evaluate", "rhs",0);

	/* cost coefficients */
	if ( !(cost = (dVector) arr_alloc(prob[1]->num->rvdOmCnt+1, double)) )
		errMsg("allocation", "evaluate", "cost", 0);
	if ( !(objxIdx = (iVector) arr_alloc(prob[1]->num->cols+1, int)) )
		errMsg("allocation", "evaluate", "objxIdx", 0);
	costTemp = expandVector(prob[1]->dBar->val, prob[1]->dBar->col, prob[1]->dBar->cnt, prob[1]->num->cols);
	for (m = 1; m <= prob[1]->num->rvdOmCnt; m++ ) {
		objxIdx[m] = prob[1]->coord->rvdOmCols[m] - 1;
		cost[m] = costTemp[prob[1]->coord->rvdOmCols[m]];
	}
	mem_free(costTemp);

	/* change the right hand side with the solution */
	chgRHSwSoln(prob[1]->bBar, prob[1]->Cbar, rhs, Xvect);

	while (3.92 * stdev > config.EVAL_ERROR * DBL_ABS(mean) || cnt < config.EVAL_MIN_ITER ) {
		/* use the stoc file to generate observations */
		generateOmega(stoc, observ+1, config.TOLERANCE, &config.EVAL_SEED[0], NULL);

		for ( m = 0; m < stoc->numOmega; m++ )
			observ[m+1] -= stoc->mean[m];          /* store the mean rv in observ */

		/* Change right-hand side with random observation */
		if ( chgRHSwObserv(subprob->lp, prob[1]->num, prob[1]->coord, observ, rhs, Xvect) ) {
			errMsg("algorithm", "evaluate", "failed to change right-hand side with random observations",0);
			return 1;
		}

		/* Change cost coefficients with random observations */
		if ( prob[1]->num->rvdOmCnt > 0 ) {
			if ( chgObjxwObserv(subprob->lp, prob[1]->num, prob[1]->coord, cost, objxIdx, observ) ) {
				errMsg("algorithm", "evaluate","failed to change cost coefficients with random observations", 0);
				return 1;
			}
		}

		changeLPSolverType(ALG_AUTOMATIC);
		if ( solveProblem(subprob->lp, subprob->name, subprob->type, &status) ) {
			if ( status == STAT_INFEASIBLE ) {
				/* subproblem is infeasible */
				printf("Warning:: Subproblem is infeasible: need to create feasibility cut.\n");
				return 1;
			}
			else {
				errMsg("algorithm", "evaluate", "failed to solve subproblem in solver", 0);
				return 1;
			}
		}

		/* use subproblem objective and compute evaluation statistics */
		obj = getObjective(subprob->lp, PROB_LP);

		if ( cnt == 0 )
			mean = obj;
		else {
			temp = mean;
			mean = mean + (obj - mean) / (double) (cnt + 1);
			variance  = (1 - 1 / (double) cnt) * variance + (cnt + 1) * (mean - temp) * (mean - temp);
			stdev = sqrt(variance/ (double) cnt);
		}
		cnt++;

		/* Print the results every once in a while for long runs */
		if (!(cnt % 100)) {
			printf(".");
			fflush(stdout);
		}
		if (!(cnt % 10000))
			printf("\nObs:%d mean:%lf   error: %lf \n0.90 CI: [%lf , %lf]\n", cnt, mean, 3.92 * stdev / mean,  mean - 1.96 * stdev, mean + 1.96 * stdev);
	}//END while loop
	mean += vXvSparse(Xvect, prob[0]->dBar);

	printEvaluationSummary(stdout, mean, stdev, cnt);
	writeEvaluationStatistics(soln, mean, stdev, cnt);

	mem_free(observ); mem_free(rhs);  mem_free(objxIdx); mem_free(cost);
	return 0;
}//END evaluate()
