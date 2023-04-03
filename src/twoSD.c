/*
 * twoSD.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */


#include <twoSD.h>

long long	MEM_USED = 0;	/* Amount of memory allocated each iteration */
cString   	outputDir;		/* output directory */
configType	config;			/* algorithm tuning parameters */

int main (int argc, char *argv[]) {
	int 	status;
	cString inputDir = NULL, probName = NULL;
	oneProblem *orig = NULL;
	timeType *tim = NULL;
	stocType *stoc = NULL;

	/* read algorithm configuration files */
	if (readConfig("../", inputDir) ) {
		errMsg("read", "main", "failed to read algorithm configuration file", 0);
		goto TERMINATE;
	}

	/* read problem information */
	parseCmdLine(argc, argv, &probName, &inputDir);

	/* open solver environment */
	openSolver();

	/* read problem SMPS input files */
	status = readFiles(inputDir, probName, &orig, &tim, &stoc);
	if ( status ) {
		errMsg("read", "main", "failed to read problem files using SMPS reader", 0);
		goto TERMINATE;
	}

	/* set up output directory: using the outputDir in config file and the input problem name */
	createOutputDir(outputDir, "twoSD", probName);

	/* launch the algorithm */
	status = algo(orig, tim, stoc, inputDir, probName);
	if ( status ) {
		errMsg("allocation", "main", "failed to solve the problem using 2-SD algorithm", 0);
		goto TERMINATE;
	}

	/* release structures and close solver environment */
	TERMINATE:
	freeOneProblem(orig);
	freeTimeType(tim);
	freeStocType(stoc);
	mem_free(probName);
	mem_free(inputDir);
	closeSolver();

	return 0;
}//END main()

void parseCmdLine(int argc, char *argv[], cString *probName, cString *inputDir) {

	for(int i=1; (i < argc); i++) {
		if ( argv[i][0] == '-' ) {
			switch ((argv[i])[1]) {
			case '?': printHelpMenu(); exit(0);
			case 'p': {
				(*probName) = (cString) arr_alloc(2*BLOCKSIZE, char);
				strcpy((*probName), argv[++i]); break;
			}
			case 'i': {
				(*inputDir) = (cString) arr_alloc(2*BLOCKSIZE, char);
				strcpy((*inputDir), argv[++i]); break;
			}
			case 'o': {
				outputDir = (cString) arr_alloc(2*BLOCKSIZE, char);
				strcpy(outputDir, argv[++i]); break;
			}
			case 'e': {
				config.EVAL_FLAG = atoi(argv[++i]);
				break;
			}
			case 'd': {
				config.DUAL_STABILITY = atoi(argv[++i]);
				break;
			}
			case 't': {
				switch (argv[++i][0]) {
				case 'l': config.EPSILON = 0.01; config.SCAN_LEN = 128; break;
				case 'n': config.EPSILON = 0.001; config.SCAN_LEN = 256; break;
				case 't': config.EPSILON = 0.0001; config.SCAN_LEN = 512; break;
				default: {
					goto TERMINATE;
					break;
				}}
				break;
			}
			case 'm': {
				config.MULTIPLE_REP = atoi(argv[++i]);
				break;
			}
			case 'c': {
				config.COMPROMISE_PROB = atoi(argv[++i]);
				break;
			}
			}
		}
		else {
			printf("Input options must begin with a '-'. Use '-?' for help.\n"); exit(0);
		}
	}

	if ( probName == NULL || inputDir == NULL || outputDir == NULL ) {
		printf("Problem name, input and output directory are mandatory input.\n");
		TERMINATE:
		if ( (*probName) ) mem_free((*probName));
		if ( outputDir ) mem_free(outputDir);
		if ( (*inputDir) ) mem_free((*inputDir));
		closeSolver(); exit(0);
	}

}//END parseCmdLine()

/* We allow only a few of the parameters to be selected through command line input. */
void printHelpMenu() {

	printf("Input options:\n");
	/* Problem name, input and output directory */
	printf("         -p string  -> problem name.\n");
	printf("         -i string  -> input directory where the problem SMPS files are saved.\n");
	printf("         -o string  -> output directory where the result files will be written.\n");
	/* EVAL_FLAG */
	printf("         -e {0,1}   -> evaluation flag which determines if the final solution will be evaluated through out of sample evalauation.\n");
	/* DUAL_STABILITY */
	printf("         -d {0,1}   -> use the dual stability test.\n");
	/* TOLERANCE */
	printf("         -t {l,n,t} -> tolerance level to be employed.\n");
	printf("                        Suggested tolerance(EPSILON, SCAN_LEN) = 'l'oose (0.01, 128), 'n'omimal (0.001, 256) and 't'ight (0.0001, 512)\n");
	/* MULTIPLE_REP */
	printf("         -m {0,1}   -> number of replications.\n");
	/* COMPROMISE_PROB */
	printf("         -c {0,1}   -> build and solve compromise problem.\n");

}//END helpMenu()

int readConfig(cString path2config, cString inputDir) {
	FILE 	*fptr;
	char	line[2*BLOCKSIZE], comment[2*BLOCKSIZE];
	int 	status, r2 = 1, maxReps = 30;

	strcpy(line, path2config);	strcat(line, "config.sd");
	fptr = fopen(line, "r");
	if ( fptr == NULL ) {
		errMsg("read", "readConfig", "failed to open configuration file", 0);
		return 1;
	}

	config.RUN_SEED = (long long *) arr_alloc(maxReps+1, long long);
	config.EVAL_SEED = (long long *) arr_alloc(maxReps+1, long long);
	config.NUM_SEEDS = 0;

	while ((status = (fscanf(fptr, "%s", line) != EOF))) {
		if (!(strcmp(line, "RUN_SEED"))) {
			fscanf(fptr, "%lld", &config.RUN_SEED[config.NUM_SEEDS+1]);
			config.NUM_SEEDS++;
			if ( config.NUM_SEEDS > maxReps + 1 ) {
				config.RUN_SEED = (long long *) mem_realloc(config.RUN_SEED, (maxReps+1)*sizeof(long long));
				maxReps *= 2;
			}
		}
		else if (!(strcmp(line, "TOLERANCE")))
			fscanf(fptr, "%lf", &config.TOLERANCE);
		else if (!(strcmp(line, "MIN_ITER")))
			fscanf(fptr, "%d", &config.MIN_ITER);
		else if (!(strcmp(line, "MAX_ITER")))
			fscanf(fptr, "%d", &config.MAX_ITER);
		else if (!(strcmp(line, "MASTER_TYPE")))
			fscanf(fptr, "%d", &config.MASTER_TYPE);
		else if (!(strcmp(line, "CUT_MULT")))
			fscanf(fptr, "%d", &config.CUT_MULT);
		else if (!(strcmp(line, "TAU")))
			fscanf(fptr, "%d", &config.TAU);
		else if (!(strcmp(line, "MIN_QUAD_SCALAR")))
			fscanf(fptr, "%lf", &config.MIN_QUAD_SCALAR);
		else if (!(strcmp(line, "MAX_QUAD_SCALAR")))
			fscanf(fptr, "%lf", &config.MAX_QUAD_SCALAR);
		else if (!(strcmp(line, "R1")))
			fscanf(fptr, "%lf", &config.R1);
		else if (!(strcmp(line, "R2")))
			fscanf(fptr, "%lf", &config.R2);
		else if (!(strcmp(line, "R3")))
			fscanf(fptr, "%lf", &config.R3);
		else if (!(strcmp(line, "DUAL_STABILITY")))
			fscanf(fptr, "%d", &config.DUAL_STABILITY);
		else if (!(strcmp(line, "PI_EVAL_START")))
			fscanf(fptr, "%d", &config.PI_EVAL_START);
		else if (!(strcmp(line, "PI_CYCLE")))
			fscanf(fptr, "%d", &config.PI_CYCLE);
		else if (!(strcmp(line, "PERCENT_PASS")))
			fscanf(fptr, "%lf", &config.PERCENT_PASS);
		else if (!(strcmp(line, "SCAN_LEN")))
			fscanf(fptr, "%d", &config.SCAN_LEN);
		else if (!(strcmp(line, "EVAL_FLAG")))
			fscanf(fptr, "%d", &config.EVAL_FLAG);
		else if (!(strcmp(line, "EVAL_SEED"))) {
			fscanf(fptr, "%lld", &config.EVAL_SEED[r2++]);
			if ( r2 > maxReps ) {
				config.RUN_SEED = (long long *) mem_realloc(config.RUN_SEED, (2*maxReps+1)*sizeof(long long));
				maxReps *= 2;
			}
		}
		else if (!(strcmp(line, "EVAL_MIN_ITER")))
			fscanf(fptr, "%d", &config.EVAL_MIN_ITER);
		else if (!(strcmp(line, "EVAL_ERROR")))
			fscanf(fptr, "%lf", &config.EVAL_ERROR);
		else if (!(strcmp(line, "PRE_EPSILON")))
			fscanf(fptr, "%lf", &config.PRE_EPSILON);
		else if (!(strcmp(line, "EPSILON")))
			fscanf(fptr, "%lf", &config.EPSILON);
		else if (!(strcmp(line, "BOOTSTRAP_REP")))
			fscanf(fptr, "%d", &config.BOOTSTRAP_REP);
		else if (!(strcmp(line, "MULTIPLE_REP")))
			fscanf(fptr, "%d", &config.MULTIPLE_REP);
		else if (!(strcmp(line, "COMPROMISE_PROB")))
			fscanf(fptr, "%d", &config.COMPROMISE_PROB);
		else if (!strcmp(line, "//"))
			fgets(comment, 2*BLOCKSIZE, fptr);
		else {
			printf ("%s\n", line);
			errMsg("read", "readConfig", "unrecognized parameter in configuration file", 1);
		}
	}

	fclose(fptr);

	config.NUM_SEEDS = minimum(config.NUM_SEEDS, r2);
	if ( config.MULTIPLE_REP > config.NUM_SEEDS ) {
		printf("Requesting to perform more replications than the number of seeds provided.\n");
		return 1;
	}

	if ( config.MULTIPLE_REP == 1 ) {
		config.COMPROMISE_PROB = 0;
	}


	return 0;
}//END readConfig()

int readFiles(cString inputDir, cString probName, oneProblem **orig, timeType **tim, stocType **stoc) {

	/* read problem core file */
	(*orig) = readCore(inputDir, probName);
	if ( (*orig) == NULL ) {
		errMsg("read", "readFiles", "failed to read problem core file", 0);
		return 1;
	}

	/* read problem time file */
	(*tim) = readTime(inputDir, probName, (*orig));
	if ( (*tim) == NULL ) {
		errMsg("read", "readFiles", "failed to read problem time file", 0);
		return 1;
	}

	(*stoc) = readStoc(inputDir, probName, (*orig), (*tim));

#ifdef INPUT_CHECK
	writeStocType((*stoc));
#endif

	return 0;
}//END readFiles()
