/*
 * main.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"
#include "cell.h"

long long	MEM_USED = 0;	/* Amount of memory allocated each iteration */
string   	outputDir;		/* output directory */
configType 	config;

int parseCmdLine(int argc, char *argv[], string algoName, string probName, string inputDir);
int twoSD(oneProblem *orig, timeType *tim, stocType *stoc, string probName);
int benders (oneProblem *orig, timeType *tim, stocType *stoc, string probName);

int main (int argc, char *argv[]) {
	int 	status;
	char 	inputDir[2*BLOCKSIZE], probName[NAMESIZE], algoName[NAMESIZE];
	oneProblem *orig = NULL;
	timeType *tim = NULL;
	stocType *stoc = NULL;

	/* open solver environment */
	openSolver();

	/* read problem information */
	if ( parseCmdLine(argc, argv, algoName, probName, inputDir) )
		goto TERMINATE;

	/* read problem SMPS input files */
	status = readFiles(inputDir, probName, &orig, &tim, &stoc);
	if ( status ) {
		errMsg("read", "main", "failed to read problem files using SMPS reader", 0);
		goto TERMINATE;
	}

	/* set up output directory: using the outputDir in config file and the input problem name */
	createOutputDir(outputDir, algoName, probName);

	/* setup the data structures based on the algorithm selected */
	if ( !(strcmp(algoName, "benders")) ) {
		if ( benders(orig, tim, stoc, probName) ) {
			errMsg("allocation", "main", "failed to solve the problem using SDDP", 0);
			goto TERMINATE;
		}
	}
	else {
		if ( twoSD(orig, tim, stoc, probName) ) {
			errMsg("allocation", "main", "failed to solve the problem using SDDP", 0);
			goto TERMINATE;
		}
	}

	/* release structures and close solver environment */
	TERMINATE:
	freeOneProblem(orig);
	freeTimeType(tim);
	freeStocType(stoc);
	closeSolver();

	return 0;
}//END main()

int parseCmdLine(int argc, char *argv[], string algoName, string probName, string inputDir) {

	outputDir = (string) arr_alloc(NAMESIZE, char);

	/* request for problem name to be solved, the path is assumed to be provided in the configuration file */
	if ( argc < 2 ) {
		printf("Please enter the algorithm you want to use: ");
		scanf("%s", algoName);
		printf("Please enter the name of the problem: ");
		scanf("%s", probName);
		strcpy(inputDir, "../spInput/");
		printf("Using default input directory: %s\n", inputDir);
		strcpy(outputDir, "../../spOutput/");
		printf("All solution files will be written to the default output directory: %s\n", outputDir);
	}
	else if ( argc < 3 ) {
		strcpy(algoName, argv[1]);
		printf("Please enter the name of the problem: ");
		scanf("%s", probName);
		strcpy(inputDir, "../spInput/");
		printf("Using default input directory: %s\n", inputDir);
		strcpy(outputDir, "../../spOutput/");
		printf("All solution files will be written to the default output directory: %s\n", outputDir);
	}
	else if ( argc < 4 ) {
		strcpy(algoName, argv[1]);
		strcpy(probName, argv[2]);
		strcpy(inputDir, "../spInput/");
		printf("Using default input directory: %s\n", inputDir);
		strcpy(outputDir, "../../spOutput/");
		printf("All solution files will be written to the default output directory: %s\n", outputDir);
	}
	else if ( argc < 5 ) {
		strcpy(algoName, argv[1]);
		strcpy(probName, argv[2]);
		strcpy(inputDir, argv[3]);
		strcpy(outputDir, "../../spOutput/");
		printf("All solution files will be written to the default output directory: %s\n", outputDir);
	}
	else {
		strcpy(algoName, argv[1]);
		strcpy(probName, argv[2]);
		strcpy(inputDir, argv[3]);
		strcpy(outputDir, argv[4]);
	}

	/* Take care of the variants in algorithm names that can be entered */
	if ( !(strcmp(algoName, "sd")) || !(strcmp(algoName, "2sd")) || !(strcmp(algoName, "SD")) || !(strcmp(algoName, "2SD")) )
		strcpy(algoName, "2sd");
	else if ( !(strcmp(algoName, "benders")) || !(strcmp(algoName, "Benders")) )
		strcpy(algoName, "benders");
	else {
		printf("Unknown algorithm name. Currently supported algorithms: \n"
				"			1. Stochastic Decomposition (2sd) \n"
				"			2. Benders Decomposition (benders)\n");
		return 1;
	}

	return 0;
}//END parseCmdLine()

void createOutputDir(string outputDir, string algoName, string probName) {
	struct stat st;
	char buffer[2*BLOCKSIZE];

	strcat(outputDir,algoName);
	strcat(outputDir,"/");
	if ( stat(outputDir, &st) ) {
		sprintf(buffer, "mkdir %s", outputDir);
		system(buffer);
	}
	strcat(outputDir, probName);
	strcat(outputDir, "/");
	if ( stat(outputDir, &st) ) {
		sprintf(buffer, "mkdir %s", outputDir);
		system(buffer);
	}
	else {
		sprintf(buffer, "rm -r %s*", outputDir);
		system(buffer);
	}

}//END createOutputDir()
