/*
 * sddp.c
 *
 *  Created on: Mar 27, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *     Contact: harsha@smu.edu
 *
 */

#include <sddp.h>

long 		MEM_USED = 0;	/* Amount of memory allocated each iteration */
string   	outputDir;		/* output directory */
configType	config;			/* algorithm tuning parameters */

int main (int argc, char *argv[]) {
	int 	status;
	char 	inputDir[2*BLOCKSIZE], probName[NAMESIZE];
	oneProblem *orig;
	timeType *tim;
	stocType *stoc;

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
	createOutputDir(outputDir, "sddp", probName);

	/* launch the algorithm */
	status = algo(orig, stoc, tim);
	if ( status ) {
		errMsg("allocation", "main", "failed to solve the problem using Time-staged SD", 0);
		goto TERMINATE;
	}

	/* release structures and close solver environment */
	TERMINATE:
	closeSolver();
	freeOneProblem(orig);
	freeTimeType(tim);
	freeStocType(stoc);

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

	fptr = fopen("config.sddp", "r");
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
		else if (!(strcmp(line, "FORWPASS_SEED")))
			fscanf(fptr, "%lld", &config.FORWPASS_SEED);
		else if (!(strcmp(line, "MAX_ITER")))
			fscanf(fptr, "%d", &config.MAX_ITER);
		else if (!(strcmp(line, "MIN_ITER")))
			fscanf(fptr, "%d", &config.MIN_ITER);
		else if (!(strcmp(line, "TOLERANCE")))
			fscanf(fptr, "%lf", &config.TOLERANCE);
		else if (!(strcmp(line, "BACKPASS_PCT")))
			fscanf(fptr, "%lf", &config.BACKPASS_PCT);
		else if (!(strcmp(line, "BACKPASS_SEED")))
			fscanf(fptr, "%lld", &config.BACKPASS_SEED);
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
