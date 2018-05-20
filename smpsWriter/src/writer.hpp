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

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <numeric>
#include <random>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <iterator>
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
int createSSNInstance(string inputDir, string outputDir);
int createTransshipInstance(string inputDir, string outputDir);

#endif /* WRITER_HPP_ */
