/*
 * main.cpp
 *
 *  Created on: Jan 28, 2018
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "writer.hpp"

void printHelpMenu();

int main (int argc, char *argv[]) {
	string probName, inputDir, outputDir;

	/* Request for problem name */
	parseCmdLine(argc, argv, &inputDir, &probName, &outputDir);

	/* Invoke the appropriate subroutine based on problem name. */
	if ( probName == "sgpf") {
		if ( createSGPFInstance(inputDir, outputDir) )
			perror("Failed to create the instance.\n");
	}
	else if (probName == "ssn" ) {
		if ( createSSNInstance(inputDir, outputDir) )
			perror("Failed to create the SSN instance.\n");
	}
	else if ( probName == "transshipment") {
		if ( createTransshipInstance(inputDir, outputDir) )
			perror("Failed to create the transshipment instance.\n");
	}
	else if ( probName == "rao") {
		if ( createRaoInstance(inputDir, outputDir, argc, argv) )
			perror("Failed to create the transshipment instance.\n");
	}
	else {
		perror("No subroutine to support this problem.\n");
	}

	return 0;
}//END main()

void parseCmdLine(int argc, char *argv[], string *inputDir, string *probName, string *outputDir) {

	for(int i=1; (i < argc); i++) {
		if ( argv[i][0] == '-' ) {
			switch ((argv[i])[1]) {
			case '?': printHelpMenu(); exit(0);
			case 'p': (*probName).assign(argv[++i]); break;
			case 'o': (*outputDir).assign(argv[++i]); break;
			case 'i': (*inputDir).assign(argv[++i]); break;
			}
		}
		else {
			printf("Input options (%s) must begin with a '-'. Use '-?' for help.\n", argv[i]); exit(0);
		}
	}

return;
}//END parseCmdLine()

void printHelpMenu() {

	cout << "Command prompt inputs." << endl;
	cout << "-i string :: Input directory" << endl;
	cout << "-o string :: Output directory" << endl;
	cout << "-p string :: Problem name" << endl;

}//END printHelpMenu()

SMPSmodel::SMPSmodel() {

	/* Default parameters */
	numPeriods = 2;
	numStages = 2;

}//END constructor()

SMPSmodel::~SMPSmodel() {

}//END destructor()

