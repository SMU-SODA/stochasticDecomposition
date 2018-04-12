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
	else {
		perror("No subroutine to support this problem.\n");
	}

	return 0;
}//END main()

void parseCmdLine(int argc, char *argv[], string *inputDir, string *probName, string *outputDir) {

	switch (argc) {
	case 4:
		*inputDir  = argv[1];
		*probName = argv[2];
		*outputDir = argv[3];
		break;
	case 3:
		*inputDir  = argv[1];
		*probName = argv[2];
		cout << "Enter an output directory to write SMPS files: ";
		cin >> *outputDir;
		break;
	case 2:
		*inputDir = argv[1];
		cout << "Enter problem name : ";
		cin >> *probName;
		cout << "Enter an output directory to write SMPS files: ";
		cin >> *outputDir;
		break;
	case 1:
		cout << "Enter input directory for data files : ";
		cin >> *probName;
		cout << "Enter problem name : ";
		cin >> *probName;
		cout << "Enter an output directory to write SMPS files: ";
		cin >> *outputDir;
		break;
	default:
		break;
	}

	return;
}//END parseCmdLine()

SMPSmodel::SMPSmodel() {

	/* Default parameters */
	numPeriods = 2;
	numStages = 2;

}//END constructor()

SMPSmodel::~SMPSmodel() {

}//END destructor()
