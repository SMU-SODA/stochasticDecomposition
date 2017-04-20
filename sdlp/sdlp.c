/*
 * sddp.c
 *
 *  Created on: Mar 27, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *     Contact: harsha@smu.edu
 *
 */

#include <sdlp.h>

long 		MEM_USED = 0;	/* Amount of memory allocated each iteration */
string   	outputDir;		/* output directory */
configType	config;			/* algorithm tuning parameters */

int main (int argc, char *argv[]) {
	int 	status;
	char 	inputDir[2*BLOCKSIZE], probName[NAMESIZE];
	oneProblem *orig = NULL;
	timeType *tim = NULL;
	stocType *stoc = NULL;

	/* open solver environment */
	openSolver();

	/* read algorithm configuration files */
	status = readConfig(inputDir);
	if ( status ) {
		errMsg("read", "main", "failed to read algorithm configuration file", 0);
		goto TERMINATE;
	}

	/* read problem information */
	/* request for problem name to be solved, the path should be provided in the configuration file */
	if ( argc < 2 )
		parseCmdLine(probName);
	else
		strcpy(probName, argv[1]);

	/* read problem SMPS input files */
	status = readFiles(inputDir, probName, &orig, &tim, &stoc);
	if ( status ) {
		errMsg("read", "main", "failed to read problem files using SMPS reader", 0);
		goto TERMINATE;
	}

	/* set up output directory: using the outputDir in config file and the input problem name */
	createOutputDir(outputDir, "sdlp", probName);

	/* launch the algorithm */
	status = algo(probName, orig, stoc, tim);
	if ( status ) {
		errMsg("allocation", "main", "failed to solve the problem using SDDP", 0);
		goto TERMINATE;
	}

	/* release structures and close solver environment */
	TERMINATE:
	freeOneProblem(orig);
	freeTimeType(tim);
	freeStocType(stoc);
	closeSolver();

	return 0;
}//END main()

void parseCmdLine(string probName) {

	printf("Please enter the name of the problem: ");
	scanf("%s", probName);

}//END parseCmdLine

int readConfig(string inputDir) {
	FILE 	*fptr;
	char	line[2*BLOCKSIZE], comment[2*BLOCKSIZE];
	int 	status;

	fptr = fopen("config.sdlp", "r");
	if ( fptr == NULL ) {
		errMsg("read", "readConfig", "failed to open configuration file", 0);
		return 1;
	}

	outputDir = (string) mem_malloc(BLOCKSIZE*sizeof(char));

	while ((status = (fscanf(fptr, "%s", line) != EOF))) {
		if (!(strcmp(line, "INPUTDIR")))
			fscanf(fptr, "%s", inputDir);
		else if (!(strcmp(line, "OUTPUTDIR")))
			fscanf(fptr, "%s", outputDir);
		else if (!(strcmp(line, "RUN_SEED")))
			fscanf(fptr, "%lld", &config.RUN_SEED);
		else if (!(strcmp(line, "MAX_ITER")))
			fscanf(fptr, "%d", &config.MAX_ITER);
		else if (!(strcmp(line, "MIN_ITER")))
			fscanf(fptr, "%d", &config.MIN_ITER);
		else if (!(strcmp(line, "TOLERANCE")))
			fscanf(fptr, "%lf", &config.TOLERANCE);
		else if (!(strcmp(line, "QUADRATIC")))
			fscanf(fptr, "%d", &config.QUADRATIC);
		else if (!(strcmp(line, "MIN_QUAD_SCALAR")))
			fscanf(fptr, "%lf", &config.MIN_QUAD_SCALAR);
		else if (!(strcmp(line, "MAX_QUAD_SCALAR")))
			fscanf(fptr, "%lf", &config.MAX_QUAD_SCALAR);
		else if (!(strcmp(line, "POLICY")))
			fscanf(fptr, "%d", &config.POLICY);
		else if (!(strcmp(line, "EVAL_FLAG")))
			fscanf(fptr, "%d", &config.EVAL_FLAG);
		else if (!(strcmp(line, "EVAL_SEED")))
			fscanf(fptr, "%lld", &config.EVAL_SEED);
		else if (!(strcmp(line, "EVAL_ERROR")))
			fscanf(fptr, "%lf", &config.EVAL_ERROR);
		else if (!(strcmp(line, "OPT_GAP")))
			fscanf(fptr, "%lf", &config.OPT_GAP);
		else if (!strcmp(line, "//"))
			fgets(comment, 2*BLOCKSIZE, fptr);
		else {
			printf ("%s\n", line);
			errMsg("read", "readConfig", "unrecognized parameter in configuration file", 1);
		}
	}

	fclose(fptr);

	return 0;
}//END readConfig()
