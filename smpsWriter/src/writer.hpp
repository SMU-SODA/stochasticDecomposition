/*
 * writer.hpp
 *
 *  Created on: Jan 28, 2018
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#ifndef WRITER_HPP_
#define WRITER_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <algorithm>
#include <string>
#include <sys/stat.h>

#include <ilcplex/ilocplex.h>

using namespace std;

#define NAMESIZE 32
#define BLOCKSIZE 256

#define TOLERANCE 0.00001

class SMPSmodel {

public:
	SMPSmodel();
	~SMPSmodel();

	int numStages, numPeriods;
	vector<string> timCols, timRows;
	vector<vector<string>> stocRows, stocCols;
	string objName;

	string dir = "../spInput";
};

void parseCmdLine(int argc, char *argv[], string *inputDir, string *probName, string *outputDir);

int createSGPFInstance(string inputDir, string outputDir);

#endif /* WRITER_HPP_ */
