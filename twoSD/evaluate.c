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

#include "utils.h"
#include "smps.h"
#include "solver.h"
#include "prob.h"
#include "multiAgentSP.h"

extern configType config;
extern string outputDir;
extern int numAgents;

/* Haven't intergrated the multi-cut and aggregated-cut versions */
void evaluateOptStoc(stocType *stoc, probType **prob, cellType **cell) {
    vector 	observ;
    double 	obj, mean, variance, stdev, temp, CI[2], totalTime;
    int		status, cnt, stat2, i, m, sOffset = 0;
    //    clock_t	tic, toc;
    FILE	*ePtr;

    if ( !(observ = (vector) arr_alloc(stoc->numOmega + 1, double)) )
        errMsg("allocation", "evaluateOpt", "observ", 0);

    printf("\n\nEvaluating optimal solution");
    //tic = clock();


    /* open a file to record evaluation results */
    ePtr = openFile(outputDir, "eval.dat", "w");

    for (i = 1; i <numAgents; i++) {

        cnt = 0.0; mean = 0.0; variance = 0.0; stdev = INFBOUND; cnt = 0;
        chgRHSwMean(prob[i]->bBar, prob[i]->Cbar, cell[i]->spRHS, cell[0]->incumbX);

        //while (3.92 * stdev > config.EVAL_ERROR * DBL_ABS(mean) || cnt < config.EVAL_MIN_ITER ) {
        //while (3.29 * stdev > config.EVAL_ERROR * DBL_ABS(mean) || cnt < config.EVAL_MIN_ITER ) {
        while (cnt < config.EVAL_MIN_ITER ) {
            /* use the stoc file to generate observations */
            generateOmega(stoc, observ, &config.EVAL_SEED);

            for ( m = 0; m < stoc->numOmega; m++ )
                observ[m] -= stoc->mean[m];          /* store the mean rv in observ */

            /* setup and solve subproblem */
            status = chgRHSwRand(cell[i]->sp->lp, prob[i]->num, prob[i]->coord, observ + sOffset -1, cell[i]->spRHS, cell[0]->incumbX, cell[i]->ID);
            if ( status ) {
                errMsg("algorithm", "evaluateOpt", "failed to setup the subproblem",0);
                exit(1);
            }

#if 0
            int     status;
            char probName[NAMESIZE];
            sprintf(probName,"evaluation_%d.lp", cell[i]->ID);
            status = writeProblem(cell[i]->sp->lp, probName);
            if ( status ) {
                errMsg("write problem", "new_subprob", "failed to write subproblems problem to file",0);
                exit(1);
            }
#endif
            status = solveProblem(cell[i]->sp->lp, cell[i]->sp->name, cell[i]->sp->type, &stat2);
            if ( status ) {
                if ( stat2 == STAT_INFEASIBLE ) {
                    /* subproblem is infeasible */
                    printf("Warning:: Subproblem is infeasible: need to create feasibility cut.\n");
                    exit(1);
                }
                else {
                    errMsg("algorithm", "evaluateOpt", "failed to solve subproblem in solver", 0);
                    exit(1);
                }
            }

            /* use subproblem objective and compute evaluation statistics */
            obj = getObjective(cell[i]->sp->lp, PROB_LP);

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
                printf("\n\nAgent%d::obs:%d mean:%lf   error: %lf \n 0.90 CI: [%lf , %lf]\n", cell[i]->ID, cnt, mean, 3.29 * stdev / mean,
                       mean - 1.645 * stdev, mean + 1.645 * stdev);
        }//END while loop
        //        toc = clock();

        //        totalTime = (toc-tic)/CLOCKS_PER_SEC;

        //        mean += vXvSparse(cell[0]->candidX, prob[0]->dBar);

        //CI[0] = mean - 1.96 * stdev;
        //CI[1] = mean + 1.96 * stdev;

        CI[0] = mean - 1.645 * stdev;
        CI[1] = mean + 1.645 * stdev;

        /* Print the value of the solution to a file and the screen */
        printf("\nAgent%d_Final Estimate :: obs:%d, mean:%lf, 0.90 C.I.: [%lf , %lf] \n\n",cell[i]->ID, cnt, mean, CI[0], CI[1]);

        fprintf(ePtr, "%d\t\t%lf\t%lf\t[%lf, %lf]\ttime:%lf\t%d\n",cell[i]->ID, mean, stdev, CI[0], CI[1], totalTime, cnt);

        if (!(prob[i]->omegas == NULL))
            sOffset += prob[i]->omegas->numRV;

    }//END agents for loop

    fclose(ePtr);
    mem_free(observ);
}

void evaluateOptSim(stocType *stoc, probType **prob, cellType **cell, simType *sim) {
    vector 	observ;
    double 	obj, mean, variance=-1.0, stdev=-1.0, temp, CI[2], totalTime, AggMean = 0.0, cx = 0.0;
    int		status, cnt, stat2, i, m, sOffset = 0;
    //    clock_t	tic, toc;
    FILE	*ePtr, *sPtr;

    if ( !(observ = (vector) arr_alloc(stoc->numOmega + 1, double)) )
        errMsg("allocation", "evaluateOpt", "observ", 0);

    printf("\n\nEvaluating optimal solution\n");
    //tic = clock();

    /* open a file to record evaluation results */
    ePtr = openFile(outputDir, "eval.dat", "w");
    sPtr = openFile(outputDir, "sol.dat", "w");

    /* set number of used observations as a given number to ensure observations used by 2-SD, MASD-Multid and MASD-Agg are the same */
    //sim->currObs = config.MAX_ITERATION;
    sim->currObs = 1000;

    for (i = 1; i < numAgents; i++) {
#ifdef EVALCHECK
        fprintf(ePtr, "agent%d:\n", i);
        //fprintf(sPtr, "agent%d:\n", i);
        printf("agent%d:\n", i);
#endif

        cnt = 0.0; mean = 0.0; variance = 0.0; stdev = INFBOUND; cnt = 0;
        chgRHSwMean(prob[i]->bBar, prob[i]->Cbar, cell[i]->spRHS, cell[0]->incumbX);

        //while (3.92 * stdev > config.EVAL_ERROR * DBL_ABS(mean) || cnt < config.EVAL_MIN_ITER ) {
        //while (3.29 * stdev > config.EVAL_ERROR * DBL_ABS(mean) || cnt < config.EVAL_MIN_ITER ) {
        while (cnt < config.EVAL_MIN_ITER ) {
            /* use the simulator to generate observations */
            simulateOmega(sim, observ);

            for ( m = 0; m < stoc->numOmega; m++ )
                observ[m] -= stoc->mean[m];          /* store the mean rv in observ */

            /* setup and solve subproblem */
            status = chgRHSwRand(cell[i]->sp->lp, prob[i]->num, prob[i]->coord, observ + sOffset -1, cell[i]->spRHS, cell[0]->incumbX, cell[i]->ID);
            if ( status ) {
                errMsg("algorithm", "evaluateOpt", "failed to setup the subproblem",0);
                exit(1);
            }

            status = solveProblem(cell[i]->sp->lp, cell[i]->sp->name, cell[i]->sp->type, &stat2);
            if ( status ) {
                if ( stat2 == STAT_INFEASIBLE ) {
                    /* subproblem is infeasible */
                    printf("Warning:: Subproblem is infeasible: need to create feasibility cut.\n");
                    exit(1);
                }
                else {
                    errMsg("evaluation", "evaluateOpt", "failed to solve subproblem in solver", 0);
                    exit(1);
                }
            }

            /* use subproblem objective and compute evaluation statistics */
            obj = getObjective(cell[i]->sp->lp, PROB_LP);

            if ( cnt == 0 )
                mean = obj;
            else {
                temp = mean;
                mean = mean + (obj - mean) / (double) (cnt + 1);
                variance  = (1 - 1 / (double) cnt) * variance
                + (cnt + 1) * (mean - temp) * (mean - temp);
                stdev = sqrt(variance/ (double) cnt);
            }
#ifdef EVALCHECK
            fprintf(ePtr, "%lf\n", obj);
            //printf("cnt:%d, obj: %lf, mean: %lf, var: %lf\n", cnt, obj, mean, variance);
#if 1
            if (i == 11){
//                if (cnt == 0) {
//                    writeProblem(cell[i]->sp->lp, "eval_ag1.lp");
//                    FILE    *sPtr0;
//                    sPtr0 = openFile(outputDir, "sol0.dat", "w");
//                    printVector(cell[0]->incumbX, cell[0]->sp->mac, sPtr0);
//                    fclose(sPtr0);
//
//                }
            	status = getPrimal(cell[i]->sp->lp, cell[i]->candidX, cell[i]->sp->mac);
            	if ( status ) {
            		errMsg("evaluation", "solveMaster", "failed to obtain the primal solution for master", 0);
            		exit(1);
            	}
            	printVector(cell[i]->candidX, cell[i]->sp->mac, sPtr);
            }
#endif


#endif
            cnt++;
            /* Print the results every once in a while for long runs */
            if (!(cnt % 100)) {
                printf(".");
                fflush(stdout);
            }
            if (!(cnt % 10000))
                printf("\n\nAgent%d::obs:%d mean:%lf   error: %lf \n 0.90 CI: [%lf , %lf]\n", cell[i]->ID, cnt, mean, 3.29 * stdev / mean, mean - 1.645 * stdev, mean + 1.645 * stdev);
                //printf("\n\nAgent%d::obs:%d mean:%lf   error: %lf \n 0.95 CI: [%lf , %lf]\n", cell[i]->ID, cnt, mean, 3.92 * stdev / mean,
                //       mean - 1.96 * stdev, mean + 1.96 * stdev);
        }//END while loop
        //        toc = clock();
        //        totalTime = (toc-tic)/CLOCKS_PER_SEC;

        //        CI[0] = mean - 1.96 * stdev;
        //        CI[1] = mean + 1.96 * stdev;

        CI[0] = mean - 1.645 * stdev;
        CI[1] = mean + 1.645 * stdev;

        /* Print the value of the solution to a file and the screen */
        printf("\nAgent%d_Final Estimate :: obs:%d, mean:%lf, 0.90 C.I.: [%lf , %lf] \n\n",cell[i]->ID, cnt, mean, CI[0], CI[1]);
        fprintf(ePtr, "%d\t\t%lf\t%lf\t[%lf, %lf]\t%d\n",cell[i]->ID, mean, stdev, CI[0], CI[1], cnt);

        if (!(prob[i]->omegas == NULL))
            sOffset += prob[i]->omegas->numRV;

        AggMean += cell[i]->weight * mean;
    }//END agents for loop

    if (numAgents < 3) {
        AggMean += vXvSparse(cell[0]->incumbX, prob[0]->dBar);
        //        CI[0] = AggMean - 1.96 * stdev;
        //        CI[1] = AggMean + 1.96 * stdev;
        CI[0] = AggMean - 1.645 * stdev;
        CI[1] = AggMean + 1.645 * stdev;
        printf("\nAgent0_Final Estimate :: mean:%lf, stdev:%lf, 0.90 C.I.: [%lf , %lf] \n\n", AggMean, stdev, CI[0], CI[1]);
        fprintf(ePtr, "0\t\t %lf\t\t %lf\t\t %lf\t\t [%lf, %lf]\n", AggMean, variance, stdev, CI[0], CI[1]);
    }
    else {
        cx=  vXvSparse(cell[0]->incumbX, prob[0]->dBar);
        AggMean = AggMean + cx;
        printf("\nAgent0_Final Estimate :: mean:%lf, cx:%lf\n\n",AggMean,cx);
        fprintf(ePtr, "0\t\t%lf\t\t%lf\n", AggMean,cx);
    }

    fclose(ePtr);
    fclose(sPtr);
    mem_free(observ);
}
