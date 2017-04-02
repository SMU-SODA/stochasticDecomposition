/*
 * evaluate.c
 *
 *  Created on: Dec 24, 2015
 *      Author: Harsha Gangammanavar
 */

#include <utils.h>
#include <smps.h>
#include <prob.h>
#include <solver.h>
#include <sddp.h>

extern configType config;

BOOL optimal(probType **prob, cellType **cell, int numStages) {
	int 	t, status, stat1, count, obs;
	double 	lb, ub, cx, mean, stdev, vari, temp, gap;

	if (cell[0]->k < config.MIN_ITER)
		return FALSE;

	/* assign optimal objective function value of the root stage obtained in the last stage as the lower bound */
	lb = cell[0]->candidEst;

	/* upper bound estimate calculation */
	cx = vXvSparse(cell[0]->candidU, prob[0]->dBar);
	/* initialize evaluation statistics */
	mean = 0.0; vari = 0.0;
	stdev = 1000000.0;
	count = 0;
	while (3.92 * stdev > config.EVAL_ERROR * DBL_ABS(mean) || count < 1000) {
		ub = 0.0;
		for ( t =  1; t < numStages; t++ ) {
			/* simulate an observation */
			obs = randInteger(&config.EVAL_SEED, cell[t]->omega->cnt);

			/* change right-hand side with endogenous information */
			computeEndoRHS(prob[t]->bBar, prob[t]->Cbar, cell[t-1]->candidU, cell[t]->rhs);

			/* change right-hand side with exogenous information */
			status = computeExoRHS(cell[t]->sp->lp, prob[t]->coord, prob[t]->num, cell[t]->omega->vals[obs], cell[t-1]->candidU, cell[t]->rhs);
			if ( status ) {
				errMsg("allocation", "forwardPass", "failed to change the right-hand side with uncertainty and state information", 0);
				return 1;
			}

			/* solve the stage problem */
			status = solveProblem(cell[t]->sp->lp, cell[t]->sp->name, PROB_LP, &stat1);
			if (status) {
				errMsg("solver", "forwardPass", "failed to solve stage problem", 0);
				return 1;
			}

			/* obtain the current primal solution */
			stat1 = getPrimal(cell[t]->sp->lp, cell[t]->candidU, prob[t]->num->cols);
			if ( stat1 ) {
				errMsg("solver", "forwardPass", "failed to obtain primal solution for the stage problem", 0);
				return 1;
			}

			/* obtain upper bound */
			ub += vXvSparse(cell[t]->candidU, prob[t]->dBar);
		}

		if ( count == 0 )
			mean = ub;
		else {
			temp = mean;
			mean = mean + (ub - mean) / (double) (count + 1);
			vari = (1 - 1 / (double) count) * vari + (count + 1) * (mean - temp) * (mean - temp);
			stdev = sqrt(vari / (double) count);
		}
		count++;
	}
	gap = cx + mean + 1.96*stdev - lb;
	printf("Lower bound = %lf\t", lb);
	printf("Upper bound = %lf (mean); %lf (stdev) with C.I. [%lf , %lf]\n", (cx+mean), stdev, (cx + mean - 1.96*stdev), (cx + mean + 1.96*stdev));
	printf("Optimality gap = %lf (%lf)\n", gap, fabs((double) gap/(cx+mean)));

	if ( gap > 0 && (fabs((double) gap/(cx+mean)) < config.OPT_GAP) )
		return TRUE;

	return FALSE;
}//END evaluate()
