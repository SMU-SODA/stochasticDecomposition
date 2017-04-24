/*
 * optimal.c
 *
 *  Created on: Apr 3, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "sdlp.h"

extern configType config;

BOOL optimal(probType **prob, cellType **cell, int numStages) {

	if ( numStages == 2 && cell[0]->k > config.MIN_ITER) {
		/* apply two-stage statistical optimality conditions */
		if (cell[1]->dualStableFlag ) {
			if ( preTest(cell[0]) ) {
				if ( fullTest(prob, cell) ) {
					cell[0]->optFlag = TRUE;
					printf("<"); fflush(stdout);
					return TRUE;
				}
				else {
					printf(">"); fflush(stdout);
				}
			}
		}
		if ( cell[0]->k >= config.MAX_ITER )
			return TRUE;
		return FALSE;
	}

	/* for multistage instances, use the trivial stopping rule. */
	if ( cell[0]->k >= config.MAX_ITER )
		return TRUE;

	return FALSE;
}//END optimal()

/* Because checking optimality is an arduous task, we first do a pre-check to determine if the full test is worthwhile. This function
 * determines whether the height at the candidate is close enough to the height at the incumbent to warrant an optimality test. */
BOOL preTest(cellType *cell) {

	/* The candidate must be within some small percentage of incumbent cut */
	/* rare situation for cell[0]->candid_est < 0 and cell[0]->incumb_est > 0 */
	/* Note: cell[0]->candidEst and cell[0]->incumbEst could be 0 */
	if (cell->candidEst >= 0){
		cell->optFlag = (cell->candidEst >= (1 - config.PRE_EPSILON) * cell->incumb->est[0]);
	}
	else
		cell->optFlag = (cell->candidEst > (1 + config.PRE_EPSILON) * cell->incumb->est[0]);

	return cell->optFlag;

}//END preTest()

BOOL fullTest(probType **prob, cellType **cell) {
	cutsType *gCuts;
	intvec  cdf, observ;
	double  est, lb = 0.0;
	int 	m, cutCnt, num_pass = 0;

	/* (a) choose good cuts */
	gCuts = chooseCuts(cell[0]->cuts, cell[0]->pi, prob[0]->num->cols);

	/* (b) calculate empirical distribution of omegas */
	if ( !(observ = (intvec) arr_alloc(cell[0]->k + 1, int)) )
		errMsg("allocation", "full_test", "failed to allocate memory to observ",0);
	if (!(cdf = (intvec) arr_alloc(cell[1]->omega->cnt, int)) )
		errMsg("allocation", "fullTest", "failed to allocate memory to cdf",0);
	empiricalDistrib(cell[1]->omega, cdf);

	for (m = 0; m < config.M; m++) {
		cutCnt = 0.0;

		/* (c) sample omegas */
		sampleOmega(cdf, observ, cell[0]->k-1);

		/* (d) reform the good cuts by plugging in the omegas */
		reformCuts(cell[1]->sigma, cell[1]->delta, cell[1]->omega, prob[1]->num, prob[1]->coord, gCuts, observ, cell[0]->k-1,
				prob[0]->lb, prob[0]->num->cols);

		/* find the highest reformed cut at the incumbent solution */
		est = maxCutHeight(gCuts, cell[0]->lb, cell[0]->k, prob[1]->coord->colsC, prob[1]->num->cntCcols, cell[0]->incumb->vals[0]);

		/* (e) Solve the master with reformed "good cuts"(all previous cuts are dropped) -> LB */
		/* In QP approach, we don't include the incumb_x * c in Sm */
		if (config.QUADRATIC)
			lb = calcTempLB(prob[0], cell[0], gCuts, cutCnt);
		else {
			est += vXvSparse(cell[0]->incumb->vals[0], prob[0]->dBar);
			printf("Lower bound calculation when solved as LP is not ready. \n");
		}


#ifdef OPT_CHECK
		printf("\nFULLTEST: replication = %d, UB = %f, LB = %f", m, AggEst, LB);
		printf("\niter = %d, replication = %d, UB = %f, LB = %f, Gap = %lf", cell[0]->k, m, AggEst, LB, DBL_ABS((AggEst - LB) / cell[0]->incumbEst));
#endif

		/* (f) compare the normalized difference btw AggEst and LB */
		/* since we know the problem is a QP problem, we don't need add the constant term c^T x \hat{x} */
		if (DBL_ABS((est - lb) / cell[0]->incumb->est[0]) <= config.OPT_GAP)
			num_pass++;

		/* (g) check No. of fails. skip out of the loop if there's no hope of meeting the condition */
		if (m + 1 - num_pass >= (1 - config.PERCENT_PASS) * config.M) {
			if(gCuts)
				freeCutsType(gCuts);
			if ( cdf) mem_free(cdf);
			if (observ) mem_free(observ);
			return FALSE;
		}
	}//END replication loop

#ifdef OPT_CHECK
	printf("\ncell[0]->k = %d, num_pass=%d\n", cell[0]->k, num_pass);
#endif


	if(gCuts)
		freeCutsType(gCuts);
	if ( cdf) mem_free(cdf);
	if (observ) mem_free(observ);
	return TRUE;
}//END fullTest()

/* This function selects all cuts whose height at the incumbent x is close to the height of the incumbent cut. These cuts together
 * are likely to provide good approximations of the recourse function at incumb_x, when they are reformed with new observations.
 * The function returns a new cut structure which contains room for cuts to be reformed. Only the _istar_ and _cut_obs_ fields of
 * each cut have been initialized. */
cutsType *chooseCuts(cutsType *cuts, vector pi, int lenX) {
	cutsType *rCuts;
	int cnt;

	rCuts = newCuts(cuts->maxCuts);

	for ( cnt = 0; cnt < cuts->cnt; cnt++ ) {
		if (pi[cuts->vals[cnt]->rowNum + 1] > 0.00001) {
			rCuts->vals[rCuts->cnt] = newCut(cuts->vals[cnt]->numIstar, cuts->vals[cnt]->numObs, lenX);
			copyIntvec(cuts->vals[cnt]->iStar, rCuts->vals[rCuts->cnt]->iStar, cuts->vals[cnt]->numIstar);
			rCuts->vals[rCuts->cnt]->rowNum = cuts->vals[cnt]->rowNum;
			rCuts->cnt++;
		}
	}

	if (rCuts->cnt == 0) {
		freeCutsType(rCuts);
		rCuts = NULL;
	}

	return rCuts;
}//END chooseCuts()

/* This function will calculate a new set of cuts based on the observations of omega passed in as _observ_, and the istar's which have already been stored
 * in the _istar_ field of each cut. If an istar field does not exist for a given observation, then a value of zero is averaged into the calculation of
 * alpha & beta. */
void reformCuts(sigmaType *sigma, deltaType *delta, omegaType *omega, numType *num, coordType *coord, cutsType *gCuts, intvec observ, int k,
		int lb, int lenX) {
	int cnt, obs, idx, count;
	iType iStar;

	/* Loop through all the cuts and reform them */
	for (cnt = 0; cnt < gCuts->cnt; cnt++) {
		/* Begin with cut coefficients of zero */
		for (idx = 0; idx <= lenX; idx++)
			gCuts->vals[cnt]->beta[idx] = 0.0;
		gCuts->vals[cnt]->alpha = 0.0;

		count = 0;
		/* Reform this cut based on resampled observations */
		for (obs = 0; obs < k; obs++) {
			/* Only sum values if the cut has an istar for this observation */
			if (observ[obs] < gCuts->vals[cnt]->numIstar) {
				iStar.sigma = gCuts->vals[cnt]->iStar[observ[obs]];
				iStar.delta = sigma->lambdaIdx[iStar.sigma];

				gCuts->vals[cnt]->alpha += sigma->vals[iStar.sigma].pib + delta->vals[iStar.delta][observ[obs]].pib;

				for (idx = 1; idx <= num->cntCcols; idx++)
					gCuts->vals[cnt]->beta[coord->colsC[idx]] += sigma->vals[iStar.sigma].piC[idx];

				for (idx = 1; idx <= num->rvColCnt; idx++)
					gCuts->vals[cnt]->beta[coord->rvCols[idx]] += delta->vals[iStar.delta][observ[obs]].piC[idx];

				count++;
			}
		}

		/* Take the average of the alpha and beta values */
		for (idx = 0; idx <= lenX; idx++)
			gCuts->vals[cnt]->beta[idx] /= (double) k;

		gCuts->vals[cnt]->alpha = (gCuts->vals[cnt]->alpha/(double) k) + (1 - (double) count / (double) k) * lb;
	}
}//END reform_cuts

/* This function forms an empirical distribution on the observations stored in omega, and calculates an integer cdf to represent the distribution.
 * An observation which has been seen n times will have n times the probability of being chosen as an observation seen only once. */
void empiricalDistrib(omegaType *omega, intvec cdf) {
	int cnt;

	/* Calculate an integer cdf distribution for observations. If the cnt is not a valid omega idx, we know that weight[cnt] is 0 */
	cdf[0] = omega->weights[0];
	for (cnt = 1; cnt < omega->cnt; cnt++)
		cdf[cnt] = cdf[cnt - 1] + omega->weights[cnt];

}//END empirical_distrib

/* This function randomly selects a new set of observations from the old set of observations stored in omega.  Entries in omega which have been observed
 * multiple times have a proportionally higher chance of being selected for the new set.  The function fills an array, assumed to be of a size equal to
 * the number of iterations which have passed, with the new set of observations. */
void sampleOmega(int *cdf, int *observ, int k) {
	int cnt, obs;
	int sample;

	/* Choose k observations according to cdf (k = number of iterations) */
	for (obs = 0; obs < k; obs++) {
		sample = randInteger(&config.EVAL_SEED, k);
		for (cnt = 0; sample > cdf[cnt]; cnt++)
			/* Loop until sample falls below cdf */;
		observ[obs] = cnt;
	}

}//END sampleOmega()

/* This function is to calculate the lower bound on the optimal value which is used in stopping rule in full_test() in optimal.c in the case of
 *  regularized approach. min c\Top x + eta + QP     s.t. Ax <= b */
double calcTempLB(probType *prob, cellType *cell, cutsType *cuts, int cutCnt) {
	double *bk; 			/* vector: b - A*incumb_x. */
	double *lambda; 		/* vector: the dual of the primal constraints. */
	double bk_lambda; 		/* scalar: bk*lambda. */
	sparseMatrix *A_Trans; 	/* sparse_matrix: the transpose of A(we call it Dbar in our code)*/
	double *A_Trans_lambda; /* vector: - A_Trans * lambda. */
	double theta; 			/* the dual of the reformed cut constraints. */
	double *Vk; 			/* the vector of scalars of cut constraints. Vk = alpha - beta * incumb_x  */
	double Vk_theta; 		/* Scalar: Vk*theta. */
	double *Bk_theta; 		/* vector: Bk_Transpose * theta, where Bk_Transpose is the matrix of cut coefficients. */
	double *q_vec; 			/* vector: c + Bk_theta - A_Trans_lambda. */
	double q_term; 			/* scalar: q_vec * q_vec. */
	double Lm; 				/* The calculated lower bound of the optimal value. */
	int cnt, i;

	if (!(bk = arr_alloc(prob->num->rows+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to bk", 0);
	if (!(lambda = arr_alloc(prob->num->rows+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to lambda", 0);
	if (!(A_Trans = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_Trans", 0);
	if (!(A_Trans->val = arr_alloc(prob->Dbar->cnt+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_Trans->val", 0);
	if (!(A_Trans->row = arr_alloc(prob->Dbar->cnt+1, int)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_Trans->row", 0);
	if (!(A_Trans->col = arr_alloc(prob->Dbar->cnt+1, int)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_Trans->col", 0);
	if (!(A_Trans_lambda = arr_alloc(prob->num->cols+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_lambda", 0);
	if (!(Vk = arr_alloc(cutCnt+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to Vk", 0);
	if (!(Bk_theta = arr_alloc(prob->num->cols+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to Bk_theta", 0);
	if (!(q_vec = arr_alloc(prob->num->cols+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to q_vec", 0);

	/* 1a. Calculate bk, which is A*incumb_x - b. Note: in fact, we are
	 ** calculating -bk here, due to the way function MSparsexvSub works. Also be aware of the one-norm. */
	for (cnt = 0; cnt < prob->num->rows; cnt++)
		bk[cnt + 1] = prob->sp->rhsx[cnt];

	/* 1b. Calculate bk = b - A * incumb_x. */
	MSparsexvSub(prob->Dbar, cell->incumb->vals[0], bk);

	/* 1c. Obtain lambda from cell->pi of original master problem constraints. */
	/* Dual values' sign need to be flipped here before assigning to lambda */
	for (cnt = 0; cnt < prob->num->rows; cnt++)
		lambda[cnt + 1] = -cell->pi[cnt + 1];

	/* 1d. Calculate bk_lambda = bk * lambda. */
	bk_lambda = vXv(bk, lambda, NULL, prob->num->rows);

	/* 2a. Calculate A_Trans */
	A_Trans->cnt = prob->Dbar->cnt;

	for (cnt = 1; cnt <= A_Trans->cnt; cnt++) {
		A_Trans->val[cnt] = prob->Dbar->val[cnt];
		A_Trans->row[cnt] = prob->Dbar->col[cnt];
		A_Trans->col[cnt] = prob->Dbar->row[cnt];
	}

	/* 2b. Calculate - A_Trans * lambda. */
	MSparsexvSub(A_Trans, lambda, A_Trans_lambda);

	/* 2c. Calculate -A_trans * lambda - c */
	for (i = 0; i < prob->num->cols; i++)
		A_Trans_lambda[i + 1] += cell->dj[i + 1];

	Vk_theta = 0.0;
	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		/* 3a. Obtain theta from c->pi */
		//MARK:cuts[1]->val[cnt]->cutObs?
		theta = ((double) (cell->k-1) / (double) cuts->vals[cnt]->numObs) * cell->pi[cuts->vals[cnt]->rowNum + 1];

		/* 3b. Obtain theta from c->pi */
		Vk_theta += theta * (cuts->vals[cnt]->alpha - vXv(cuts->vals[cnt]->beta, cell->incumb->vals[0], NULL, prob->num->cols));

		/* 3c. Calculate Bk_theta = theta * beta */
		for (i = 1; i <= prob->num->cols; i++)
			Bk_theta[i] += theta * cuts->vals[cnt]->beta[i];
	}

	/* 4. Calculate the quadratic vector q_vec = Bk_theta[i] - A_Trans_lambda[i].*/
	for (i = 1; i <= prob->num->cols; i++)
		q_vec[i] = prob->dBar->val[i] - Bk_theta[i] - A_Trans_lambda[i];

	q_term = vXv(q_vec, q_vec, NULL, prob->num->cols);

	/* 5. Calculate the lower bound*/
	Lm = Vk_theta + bk_lambda - q_term / cell->incumb->quadScalar/ 2.0;

	mem_free(bk);
	mem_free(lambda);
	mem_free(A_Trans->col);
	mem_free(A_Trans->row);
	mem_free(A_Trans->val);
	mem_free(A_Trans);
	mem_free(A_Trans_lambda);
	mem_free(Vk);
	mem_free(Bk_theta);
	mem_free(q_vec);

	return Lm;
}//END calcTempLB()
