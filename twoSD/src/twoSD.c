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
void parseCmdLine(int argc, char *argv[], cString *probName, cString *inputDir);
void printHelpMenu();

long long	MEM_USED = 0;	/* amount of memory allocated each iteration */
cString	outputDir;			/* output directory */
configType	config;			/* algorithm tuning parameters */

int main (int argc, char *argv[]) {
	int 	status;
	cString inputDir = NULL, probName = NULL;
	oneProblem *orig = NULL;
	timeType *tim = NULL;
	stocType *stoc = NULL;
	
	outputDir = NULL;
	/* read the default algorithm configuration parameters */
	if (readConfig("C:\\Users\\stabr\\Documents\\GitHub\\stochasticDecomposition\\twoSD\\src\\", inputDir)) {
//	if (readConfig("./src/", inputDir)) {
		errMsg("read", "main", "failed to read algorithm configuration file", 0);
		goto TERMINATE;
	}


	/* read problem information */
	parseCmdLine(argc, argv, &probName, &inputDir);
	
	/* read problem SMPS input files */
	status = readFiles(inputDir, probName, &orig, &tim, &stoc);
	if ( status ) {
		errMsg("read", "main", "failed to read problem files using SMPS reader", 0);
		goto TERMINATE;
	}
	
	/* set up output directory: using the outputDir in config file and the input problem name */
	createOutputDir(outputDir, "twoSD-integer", probName);
	
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
	closeSolver();
	mem_free(probName); mem_free(inputDir);

	return 0;
}//END main()

void parseCmdLine(int argc, char *argv[], cString *probName, cString *inputDir) {


	for(int i=1; (i < argc); i++) {
		if ( argv[i][0] == '-' ) {
			switch ((argv[i])[1]) {
			case '?': printHelpMenu(); exit(0);
			case 'r': { //getting the probName form the user
				(*probName) = (cString)arr_alloc(2 * BLOCKSIZE, char);
				printf("Please enter the problem name: ");
				scanf("%s", (*probName)); break;
			}
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
			case 's':{
				config.SAMPLE_INCREMENT = atoi(argv[++i]);
				break;
			}
			}
		}
		else {
			printf("Input options must begin with a '-'. Use '-?' for help.\n"); exit(0);
		}
	}


	if (probName == NULL || inputDir == NULL || outputDir == NULL) {
		printf("Problem name, input and output directory are mandatory input.\n");
	TERMINATE:
		if ((*probName)) mem_free((*probName));
		if (outputDir) mem_free(outputDir);
		if ((*inputDir)) mem_free((*inputDir));
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
	printf("         -m {0,1}   -> use multiple replication.\n");
	/* COMPROMISE_PROB */
	printf("         -c {0,1}   -> build and solve compromise problem.\n");

}//END helpMenu()
