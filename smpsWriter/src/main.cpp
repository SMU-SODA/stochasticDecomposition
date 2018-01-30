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
	string probName;

	/* Request for problem name */
	parseCmdLine(argc, argv, &probName);

	/* Invoke the appropriate subroutine based on problem name. */
	if ( probName == "sgpf") {
		if ( createSGPFInstance() )
			perror("Failed to create the instance.\n");
	}
	else {
		perror("No subroutine to support this problem.\n");
	}

	return 0;
}//END main()

void parseCmdLine(int argc, char *argv[], string *probName) {

	switch (argc) {
	case 2:
		*probName = argv[1];
		break;
	case 1:
		cout << "Enter problem name : ";
		cin >> *probName;
		break;
	default:
		break;
	}

	return;
}//END parseCmdLine()
