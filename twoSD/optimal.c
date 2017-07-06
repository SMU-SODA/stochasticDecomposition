/*
 * optimal.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

//
//  optimal.c
//  multiAgentSP
//
//  Created by Shasha Wang on 3/28/16.
//  Copyright © 2016 Shasha Wang. All rights reserved.
//

#include "utils.h"
#include "smps.h"
#include "solver.h"
#include "prob.h"
#include "multiAgentSP.h"

extern configType config;
extern int numAgents;

/* This function determines whether or not the current incumbent solution is considered to be optimal. Optimality is guarenteed if the
 * following criteria are satisfied:
 * 		0. Minimum number of iterations have been completed.
 * 		1. Dual solution set has stabilized for an agent.
 * 		2. If dual solution set is stable, the pre-test checks for "convergence" of objective function estimates at each agent.
 * 			- If estimates are satisfactory, then the agent is switched off. No more problems solved for that agent.
 * 		3. Full test is based on boot-strapping, and checks the gap between primal (upper) and dual (lower) values.
 * The pre-test is performed only after the dual solution set has stabilized, and the full test is performed only if all the agents have
 * satisfied the pre-test. */
BOOL optimal(probType **prob, cellType **cell){
	int n;

	/* ensure that the minimum number of iterations have been completed */
	if (cell[0]->k > config.MIN_ITERATION) {
        cell[0]->optFlag = TRUE;
		for ( n = 1; n < numAgents; n++ ) {
			/* check to see if the dual solution set is stable */
			cell[n]->optFlag = FALSE;
            if ( cell[n]->dualStableFlag ) {
				/* perform the pre-test */
				preTest(cell[n]);
            }
			cell[0]->optFlag &= cell[n]->optFlag;
		}

		/* If all agents have satisfied pre-test then we perform the full test */
		if ( cell[0]->optFlag ) {
            printf ("\n(>%d)", cell[0]->k);
			if (full_test(prob, cell)) {
				/* full test satisfied */
				cell[0]->optFlag = TRUE;
				printf ("\n(*%d)", cell[0]->k);
				return TRUE;
			}
            else{
                cell[0]->optFlag = FALSE;
                printf ("\n(<<%d)", cell[0]->k);
            }
		}
	}
	return FALSE;
}//optimal()

/* Because checking optimality is an arduous task, we first do a pre-check to determine if the full test is worthwhile. This function
 * determines whether the height at the candidate is close enough to the height at the incumbent to warrant an optimality test. */
BOOL preTest(cellType *cell) {

	/* The candidate must be within some small percentage of incumbent cut */
	/* rare situation for cell[0]->candid_est < 0 and cell[0]->incumb_est > 0 */
	/* Note: cell[0]->candidEst and cell[0]->incumbEst could be 0 */
	if (cell->candidEst >= 0){
		cell->optFlag = (cell->candidEst >= (1 - config.PRE_EPSILON) * cell->incumbEst);
	}
	else
		cell->optFlag = (cell->candidEst > (1 + config.PRE_EPSILON) * cell->incumbEst);

	return cell->optFlag;

}//END preTest()

/* This function performs a complete statistical test of optimality. First, it selects cuts whose height at the incumbent is "close" to
 * incumbent cut's height.  Then, it performs M resamplings of the observations in omega, and reforms the selected cuts with respect
 * to these observations (as if *they* were observed instead of the actual omega).  For each of the M resamplings, a master program
 * containing the reformed cuts is solved, and if almost all of the solutions to these master programs.
 * Multi-cut version selects "good" cuts within each agent, and reforms theses cuts. Then the cut value are aggregated.
 * Single-cut version first selects "good" cuts in the master prob. Then by using the information of these "good" cuts to find the
 * corresponding cuts in each agent. After reforming these agent cuts, aggregated these reformed cuts as one single cut and add it to the
 * master problem */
BOOL full_test(probType **prob, cellType **cell){
	cutType **gCuts;
	vector  Est;
	double  ht, AggEst, LB=0.0, error_sum = 0.0;
	intvec  *cdf;
	int *observ, m, i, j, cutCnt, num_pass = 0, num_failed = 0;

	if ( !(gCuts = (cutType **) arr_alloc(numAgents, cutType *)) )
		errMsg("allocation", "full_test", "failed to allocate memory to gCuts",0);
	gCuts[0] = NULL;

	if (!(cdf = (intvec *) arr_alloc(numAgents, intvec)) )
		errMsg("allocation", "full_test", "failed to allocate memory to cdf",0);

	for (i = 1; i < numAgents; i++) {
		/* (a) choose good cuts */
		gCuts[i] = chooseCuts(prob[0], cell[i], cell[0]);

		/* (b) calculate empirical distribution of omegas */
		if ( !(cdf[i] = (intvec) arr_alloc(cell[i]->omega->cnt+1, int)) )
			errMsg("allocation", "full_test", "failed to allocate memory to cdf",0);

		empirical_distrib(cell[i]->omega, cdf[i]);
	}

	if ( config.MULTI_CUT ) {
		/* holds estimates for each agent to be used later for aggregation, so need one for each agent */
		if ( !(Est = (vector) arr_alloc(numAgents, double)) )
			errMsg("allocation", "full_test", "failed to allocate memory to Est",0);
	}
	else {
		/* holds aggregated estimate from each cut to be later used for maximization, so need one for each "good cut". */
		if ( !(Est = (vector) arr_alloc(gCuts[1]->cnt, double)) )
			errMsg("allocation", "full_test", "failed to allocate memory to Est",0);
	}

	for (m = 0; m < config.M; m++) {
		AggEst = 0.0;
		cutCnt = 0.0;

        if (!(config.MULTI_CUT)){
            for (j = 0; j < gCuts[1]->cnt; j++)
                Est[j] = 0.0;
        }

		for (i = 1; i < numAgents; i++) {
			if (!(gCuts[i] == NULL)){
				if ( !(observ = (intvec) arr_alloc(cell[0]->k + 1, int)) )
					errMsg("allocation", "full_test", "failed to allocate memory to observ",0);

				/* (c) sample omegas */
				sampleOmega(cdf[i], observ, cell[0]->k-1);

				/* (d) reform the good cuts by plugging in the
                 omegas */
				reform_cuts(cell[i]->sigma, cell[i]->delta, cell[i]->omega, prob[i]->num, prob[i]->coord, gCuts[i], observ, cell[0]->k-1, cell[0]->lbType, prob[0]->lb, prob[0]->num->cols);

				if ( config.MULTI_CUT ) {
					/* (e1) find out the best reformed cut -> Est at the incumbent solution */
					Est[i] = gCuts[i]->val[0]->alpha - vXv(gCuts[i]->val[0]->beta, cell[0]->incumbX, NULL, prob[0]->num->cols);
					for (j = 1; j < gCuts[i]->cnt; j++) {
						ht = gCuts[i]->val[j]->alpha - vXv(gCuts[i]->val[j]->beta, cell[0]->incumbX, NULL, prob[0]->num->cols);
						if (Est[i] < ht)
							Est[i] = ht;
					}
					/* (e2) aggregate the estimate of each agent */
					AggEst += Est[i] * cell[i]->weight;
					cutCnt += gCuts[i]->cnt;        /* cutCnt is used to allocate memory in cal_temp_lb() */
				}
				else {
					for (j = 0; j < gCuts[i]->cnt; j++)
						Est[j] += cell[i]->weight* (gCuts[i]->val[j]->alpha - vXv(gCuts[i]->val[j]->beta, cell[0]->incumbX, NULL, prob[0]->num->cols));
				}

				mem_free(observ);
			}
		}//END agent loop

        /* for single cuts version */
		if ( !(config.MULTI_CUT) ) {
			AggEst = Est[0];
            //MARK: why use Est[i] = 0 here?
            //Est[i] = 0;
			for ( j = 1; j < gCuts[1]->cnt; j++ )
				if ( Est[j] > AggEst )
					AggEst = Est[j];
			cutCnt = gCuts[1]->cnt;
		}

		/* (f) Solve the master with reformed "good cuts"(all previous cuts are dropped) -> LB */
		/* In QP approach, we don't include the incumb_x * c in Sm */
		if (config.MASTERTYPE == 1)
			AggEst += vXvSparse(cell[0]->incumbX, prob[0]->dBar);
		else
			LB = cal_temp_lb(prob[0], cell, gCuts, cutCnt);


#ifdef OPT_CHECK
		printf("\nFULLTEST: replication = %d, UB = %f, LB = %f", m, AggEst, LB);
#endif

		printf("\niter = %d, replication = %d, UB = %f, LB = %f, Gap = %lf", cell[0]->k, m, AggEst, LB, DBL_ABS((AggEst - LB) / cell[0]->incumbEst));
		/* (g) compare the normalized difference btw AggEst and LB */
		/* since we know the problem is a QP problem, we don't need add the constant term c^T x \hat{x} */
		if (DBL_ABS((AggEst - LB) / cell[0]->incumbEst) <= config.OPTGAP)
			num_pass++;
		else {
			/* Sum up errors of failed replications in full test. Only at Max_iter. */
			if (cell[0]->k >= config.MAX_ITERATION) {
				num_failed++;
				error_sum += (AggEst - LB) / cell[0]->incumbEst;
			}
		}

		/* (h) if the flag for boostrap test is disabled, then we simply claimed all tests passed the test */
		if (config.BOOTSTRAP_TEST == 0)
			num_pass = config.M;

		/* (i1) check No. of fails. skip out of the loop if there's no hope of meeting the condition */
		if (m + 1 - num_pass >= (1 - config.PERCENT_PASS) * config.M) {
			/* Record average error of failed replications in full test. Only at Max_iter */
			if (cell[0]->k >= config.MAX_ITERATION)
				cell[0]->full_test_error = error_sum / num_failed;
			goto TERMINATE;
		}
	}//END replication loop

#ifdef OPT_CHECK
	printf("\ncell[0]->k = %d, num_pass=%d\n", cell[0]->k, num_pass);
#endif


	mem_free(Est);
	freeIntvec(cdf,numAgents);
	if(gCuts)
		freecutType(gCuts);
	return TRUE;

	/* continue to re-forming cuts such that the variance can be calculated */
TERMINATE:
	mem_free(Est);
	freeIntvec(cdf,numAgents);
	if(gCuts)
		freecutType(gCuts);
	return FALSE;
}//END full_test()

/* This function selects all cuts whose height at the incumbent x is close to the height of the incumbent cut. These cuts together
 * are likely to provide good approximations of the recourse function at incumb_x, when they are reformed with new observations.
 * The function returns a new cut structure which contains room for cuts to be reformed. Only the _istar_ and _cut_obs_ fields of
 * each cut have been initialized. */
cutType *chooseCuts(probType *prob, cellType *cell, cellType *master) {
	cutType *cuts;
	int cnt;

	cuts = newCuts(cell->maxCuts);

	if ( config.MULTI_CUT ) {
		for (cnt = 0; cnt < cell->cuts->cnt; cnt++) {
			/* Choosing cuts with nonzero dual multipliers (indicates these cuts are binding) */
			if (master->pi[cell->cuts->val[cnt]->rowNum + 1] > 0.00001) {
				cuts->val[cuts->cnt] = newCut(prob->num->cols, cell->cuts->val[cnt]->omegaCnt, cell->cuts->val[cnt]->cutObs);
				copyIntvec(cell->cuts->val[cnt]->iStar, cuts->val[cuts->cnt]->iStar, cell->cuts->val[cnt]->omegaCnt);
				cuts->val[cuts->cnt]->rowNum = cell->cuts->val[cnt]->rowNum;
				cuts->cnt++;
			}
		}
	}
	else {
		for ( cnt = 0; cnt < master->cuts->cnt; cnt++ ) {
			if (master->pi[master->cuts->val[cnt]->rowNum + 1] > 0.00001) {
				cuts->val[cuts->cnt] = newCut(prob->num->cols, cell->cuts->val[cnt]->omegaCnt, cell->cuts->val[cnt]->cutObs);
				copyIntvec(cell->cuts->val[cnt]->iStar, cuts->val[cuts->cnt]->iStar, cell->cuts->val[cnt]->omegaCnt);
				cuts->val[cuts->cnt]->rowNum = master->cuts->val[cnt]->rowNum;
				cuts->cnt++;
			}
		}
	}

	if (cuts->cnt == 0) {
		freeCuts(cuts);
		cuts = NULL;
	}

	return cuts;
}//END choose_cuts

/***********************************************************************\
 ** This function forms an empirical distribution on the observations
 ** stored in omega, and calculates an integer cdf to represent the
 ** distribution.  An observation which has been seen n times will have n
 ** times the probability of being chosen as an observation seen only once.
 \***********************************************************************/
void empirical_distrib(omegaType *omega, int *cdf) {
	int cnt;

	/* Calculate an integer cdf distribution for observations */
	/* If the cnt is not a valid omega idx, we know that weight[cnt] is 0 */
	cdf[0] = omega->weight[0];
	for (cnt = 1; cnt < omega->cnt; cnt++)
		cdf[cnt] = cdf[cnt - 1] + omega->weight[cnt];
}//END empirical_distrib

/***********************************************************************\
 ** This function randomly selects a new set of observations from the
 ** old set of observations stored in omega.  Entries in omega which
 ** have been observed multiple times have a proportionally higher
 ** chance of being selected for the new set.  The function fills an
 ** array, assumed to be of a size equal to the number of iterations
 ** which have passed, with the new set of observations.
 ***********************************************************************/
void sampleOmega(int *cdf, int *observ, int k) {
	int cnt, obs;
	int sample;

	/* Choose k observations according to cdf (k = number of iterations) */
	for (obs = 0; obs < k; obs++) {
		sample = randFun(k, &config.EVAL_SEED1);
		for (cnt = 0; sample > cdf[cnt]; cnt++)
			/* Loop until sample falls below cdf */;
		observ[obs] = cnt;
	}
}//END sampleOmega

/***********************************************************************\
 ** This function will calculate a new set of cuts based on the
 ** observations of omega passed in as _observ_, and the istar's
 ** which have already been stored in the _istar_ field of each cut.
 ** If an istar field does not exist for a given observation, then
 ** a value of zero is averaged into the calculation of alpha & beta.
 \***********************************************************************/
void reform_cuts(sigmaType *sigma, deltaType *delta, omegaType *omega, numType *num, coordType *coord, cutType *gCuts, int *observ, int k, int lbType, int lb, int lenX) {
	int cnt, obs, idx, count;
	iType iStar;

	/* Loop through all the cuts and reform them */
	for (cnt = 0; cnt < gCuts->cnt; cnt++) {
		/* Begin with cut coefficients of zero */
		for (idx = 0; idx <= lenX; idx++)
			gCuts->val[cnt]->beta[idx] = 0.0;
		gCuts->val[cnt]->alpha = 0.0;

		count = 0;
		/* Reform this cut based on resampled observations */
		for (obs = 0; obs < k; obs++) {
			/* Only sum values if the cut has an istar for this observation */
			if (observ[obs] < gCuts->val[cnt]->omegaCnt) {
				iStar.sigma = gCuts->val[cnt]->iStar[observ[obs]];
				iStar.delta = sigma->lambdaIdx[iStar.sigma];

				gCuts->val[cnt]->alpha += sigma->vals[iStar.sigma].b + delta->vals[iStar.delta][observ[obs]].b;

				for (idx = 1; idx <= num->cntCcols; idx++)
					gCuts->val[cnt]->beta[coord->colsC[idx]] += sigma->vals[iStar.sigma].C[idx];

				for (idx = 1; idx <= num->rvColCnt; idx++)
					gCuts->val[cnt]->beta[coord->rvCols[idx]] += delta->vals[iStar.delta][observ[obs]].C[idx];

				count++;
			}
		}

		/* Take the average of the alpha and beta values */
		for (idx = 0; idx <= lenX; idx++)
			gCuts->val[cnt]->beta[idx] /= (double) k;

		gCuts->val[cnt]->alpha /= (double) k;

		if (lbType == NONTRIVIAL)
			gCuts->val[cnt]->alpha += (1 - (double) count / (double) k) * lb;
	}
}//END reform_cuts

/**********************************************************************
 ** This function returns a uniform random number between [0, greatest-1]
 ** using our own random number generator.  This is not so good, since
 ** we are all running off the same seed... we ought to have different
 ** streams of random numbers.
 **********************************************************************/
int randFun(int greatest, long long *seed) {
	return (int) (randUniform(seed) * greatest);
}//END randFun

/****************************************************************************\
 This function is to calculate the lower bound on the optimal value which
 is used in stopping rule in full_test() in optimal.c in the case of
 regularized approach. min c\Top x + eta + QP     s.t. Ax <= b
 \****************************************************************************/
double cal_temp_lb(probType *p, cellType **c, cutType **cuts, int cutCnt) {
	double *bk; /* vector: b - A*incumb_x. */
	double *lambda; /* vector: the dual of the primal constraints. */
	double bk_lambda; /* scalar: bk*lambda. */
	sparseMatrix *A_Trans; /* sparse_matrix: the transpose of A(we call it Dbar in our code)*/
	double *A_Trans_lambda; /* vector: - A_Trans * lambda. */
	double theta; /* the dual of the reformed cut constraints. */
	double *Vk; /* the vector of scalars of cut constraints. Vk = alpha - beta * incumb_x  */
	double Vk_theta; /* Scalar: Vk*theta. */
	double *Bk_theta; /* vector: Bk_Transpose * theta, where Bk_Transpose is the matrix of cut coefficients. */
	double *q_vec; /* vector: c + Bk_theta - A_Trans_lambda. */
	double q_term; /* scalar: q_vec * q_vec. */
	double Lm; /* The calculated lower bound of the optimal value. */
	int cnt, i, n;

	if (!(bk = arr_alloc(p->num->rows+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to bk", 0);
	if (!(lambda = arr_alloc(p->num->rows+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to lambda", 0);
	if (!(A_Trans = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_Trans", 0);
	if (!(A_Trans->val = arr_alloc(p->Dbar->cnt+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_Trans->val", 0);
	if (!(A_Trans->row = arr_alloc(p->Dbar->cnt+1, int)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_Trans->row", 0);
	if (!(A_Trans->col = arr_alloc(p->Dbar->cnt+1, int)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_Trans->col", 0);
	if (!(A_Trans_lambda = arr_alloc(p->num->cols+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_lambda", 0);
	if (!(Vk = arr_alloc(cutCnt+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to Vk", 0);
	if (!(Bk_theta = arr_alloc(p->num->cols+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to Bk_theta", 0);
	if (!(q_vec = arr_alloc(p->num->cols+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to q_vec", 0);

	/* 1a. Calculate bk, which is A*incumb_x - b. Note: in fact, we are
	 ** calculating -bk here, due to the way function MSparsexvSub works. Also be aware of the one-norm. */
	for (cnt = 0; cnt < p->num->rows; cnt++)
		bk[cnt + 1] = p->sp->rhsx[cnt];

	/* 1b. Calculate bk = b - A * incumb_x. */
	MSparsexvSub(p->Dbar, c[0]->incumbX, bk);

	/* 1c. Obtain lambda from cell->pi of original master problem constraints. */
	/* Dual values' sign need to be flipped here before assigning to lambda */
	for (cnt = 0; cnt < p->num->rows; cnt++)
		lambda[cnt + 1] = -c[0]->pi[cnt + 1];

	/* 1d. Calculate bk_lambda = bk * lambda. */
	bk_lambda = vXv(bk, lambda, NULL, p->num->rows);

	/* 2a. Calculate A_Trans */
	A_Trans->cnt = p->Dbar->cnt;

	for (cnt = 1; cnt <= A_Trans->cnt; cnt++) {
		A_Trans->val[cnt] = p->Dbar->val[cnt];
		A_Trans->row[cnt] = p->Dbar->col[cnt];
		A_Trans->col[cnt] = p->Dbar->row[cnt];
	}

	/* 2b. Calculate - A_Trans * lambda. */
	MSparsexvSub(A_Trans, lambda, A_Trans_lambda);

	/* 2c. Calculate -A_trans * lambda - c */
	for (i = 0; i < p->num->cols; i++)
		A_Trans_lambda[i + 1] += c[0]->di[i + 1];

	Vk_theta = 0.0;
	if ( config.MULTI_CUT ) {
		for ( n = 1; n < numAgents; n++) {
			/* 3. Calculate Vk_theta = Vk * theta. */
			if (!(cuts[n] == NULL)){
				for (cnt = 0; cnt < cuts[n]->cnt; cnt++) {
					/* 3a. Obtain theta from c->pi */
					theta = ((double) (c[0]->k-1) / (double) cuts[n]->val[cnt]->cutObs) * c[0]->pi[cuts[n]->val[cnt]->rowNum + 1];

					Vk_theta += theta * (cuts[n]->val[cnt]->alpha - vXv(cuts[n]->val[cnt]->beta, c[0]->incumbX, NULL, p->num->cols));

					/* 3b. Calculate Bk_theta = theta * beta */
					for (i = 1; i <= p->num->cols; i++)
						Bk_theta[i] += theta * cuts[n]->val[cnt]->beta[i];
				}
			}
		}
	}
	else {
		for (cnt = 0; cnt < cuts[1]->cnt; cnt++) {
			/* 3a. Obtain theta from c->pi */
            //MARK:cuts[1]->val[cnt]->cutObs?
			theta = ((double) (c[0]->k-1) / (double) cuts[1]->val[cnt]->cutObs) * c[0]->pi[cuts[1]->val[cnt]->rowNum + 1];
			for ( n = 1; n < numAgents; n++) {
				/* 3a. Obtain theta from c->pi */
				Vk_theta += theta * c[n]->weight*(cuts[n]->val[cnt]->alpha - vXv(cuts[n]->val[cnt]->beta, c[0]->incumbX, NULL, p->num->cols));

				/* 3b. Calculate Bk_theta = theta * beta */
				for (i = 1; i <= p->num->cols; i++)
					Bk_theta[i] += theta * c[n]->weight * cuts[n]->val[cnt]->beta[i];
			}
		}
	}

	/* 4. Calculate the quadratic vector q_vec = Bk_theta[i] - A_Trans_lambda[i].*/
	for (i = 1; i <= p->num->cols; i++)
		q_vec[i] = p->dBar->val[i] - Bk_theta[i] - A_Trans_lambda[i];

	q_term = vXv(q_vec, q_vec, NULL, p->num->cols);

	/* 5. Calculate the lower bound*/
	Lm = Vk_theta + bk_lambda - q_term / c[0]->quadScalar / 2.0;

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
}//END cal_temp_lb

double calc_var1(double *x, double *mean_value, double *stdev_value, int length) {
	double mean, vari, temp;
	int count;
	double stdev;
	stdev = 10000000.0;
	temp = 0.0;
	mean = x[0];
	vari = 0.0;

	for (count = 1; count < length; count++) {
		temp = mean;
		mean = mean + (x[count] - mean) / (double) (count + 1);
		vari = (1 - 1 / (double) count) * vari
				+ (count + 1) * (mean - temp) * (mean - temp);
	}

	if (mean_value != NULL)
	{
		*mean_value = mean;
	}
	if (stdev_value != NULL)
	{
		stdev = sqrt(vari / (double) count);
		*stdev_value = stdev;
	}

	return vari;

}//END calc_var1

double cal_lmn(cellType *c, cutType *T, vector pi, vector incumbX, int lenX){
	double theta=0.0, Vk_theta=0.0;
	int cnt;

	for (cnt = 0; cnt < T->cnt; cnt++) {
		/* 3a. Obtain theta from c->pi */
		theta = ((double) c->k / (double) T->val[cnt]->cutObs) * pi[T->val[cnt]->rowNum + 1];

		Vk_theta += theta * (T->val[cnt]->alpha - vXv(T->val[cnt]->beta, incumbX, NULL, lenX));
	}

	return Vk_theta;
}

void freeIntvec(intvec *cdf,int length) {
	int i;

	for (i = 1; i < length; i++)
		mem_free(cdf[i]);
	mem_free(cdf);
}//END freeIntvec

void freeArray(double **x, int length) {
	int i;

	for (i = 1; i < length; i++)
		mem_free(x[i]);
	mem_free(x);
}//END freeArray

void freecutType(cutsType **cut) {
	int n;

	for ( n = 0; n < numAgents; n++ ) {
		if(cut[n])
			freeCuts(cut[n]);
	}

	mem_free(cut);

}//END freecutType()

