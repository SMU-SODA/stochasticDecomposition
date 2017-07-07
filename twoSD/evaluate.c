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

int evaluate(stocType *stoc, probType **prob, cellType *cell, vector Xvect) {
	vector 	observ, rhs;
	double 	obj, mean, variance, stdev, temp;
	int		cnt, status, m;

	if ( !(observ = (vector) arr_alloc(stoc->numOmega + 1, double)) )
		errMsg("allocation", "evaluateOpt", "observ", 0);

	printf("\n\nEvaluating optimal solution");

	/* initialize parameters used for evaluations */
	cnt = 0.0; mean = 0.0; variance = 0.0; stdev = INFBOUND; cnt = 0;
	if (!(rhs =(vector) arr_alloc(prob[1]->num->rows+1, double)))
		errMsg("Allocation", "computeRhs", "rhs",0);

	/* change the right hand side with the solution */
	chgRHSwSoln(prob[1]->bBar, prob[1]->Cbar, rhs, Xvect);
	while (3.92 * stdev > config.EVAL_ERROR * DBL_ABS(mean) || cnt < config.EVAL_MIN_ITER ) {
		/* use the stoc file to generate observations */
		generateOmega(stoc, observ, &config.EVAL_SEED);

		for ( m = 0; m < stoc->numOmega; m++ )
			observ[m] -= stoc->mean[m];          /* store the mean rv in observ */

		/* Change right-hand side with random observation */
		if ( chgRHSwObserv(cell->subprob->lp, prob[1]->num, prob[1]->coord, observ-1, cell->spRHS, cell->incumbX) ) {
			errMsg("algorithm", "evaluateOpt", "failed to setup the subproblem",0);
			return 1;
		}

		if ( solveProblem(cell->subprob->lp, cell->subprob->name, cell->subprob->type, &status) ) {
			if ( status == STAT_INFEASIBLE ) {
				/* subproblem is infeasible */
				printf("Warning:: Subproblem is infeasible: need to create feasibility cut.\n");
				return 1;
			}
			else {
				errMsg("algorithm", "evaluateOpt", "failed to solve subproblem in solver", 0);
				return 1;
			}
		}

		/* use subproblem objective and compute evaluation statistics */
		obj = getObjective(cell->subprob->lp, PROB_LP);

		if ( cnt == 0 )
			mean = obj;
		else {
			temp = mean;
			mean = mean + (obj - mean) / (double) (cnt + 1);
			variance  = (1 - 1 / (double) cnt) * variance
					+ (cnt + 1) * (mean - temp) * (mean - temp);
			stdev = sqrt(variance/ (double) cnt);
		}
		cnt++;
		/* Print the results every once in a while for long runs */
		if (!(cnt % 100)) {
			printf(".");
			fflush(stdout);
		}
		if (!(cnt % 10000))
			printf("\n\nObs:%d mean:%lf   error: %lf \n 0.90 CI: [%lf , %lf]\n", cnt, mean, 3.29 * stdev / mean,  mean - 1.645 * stdev, mean + 1.645 * stdev);
	}//END while loop

	mem_free(observ); mem_free(rhs);
	return 0;

}//END evaluate()
