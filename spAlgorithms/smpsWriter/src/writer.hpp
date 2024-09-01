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

#define CPP_TOL	0.00001
#define NAMESIZE 32
#define BLOCKSIZE 256

using namespace std;

class SMPSmodel {

public:
	SMPSmodel();
	SMPSmodel(int t);
	~SMPSmodel();

	int numStages, numPeriods;
	vector<string> timeCols, timeRows;
	vector<vector<string>> stocRows, stocCols;
	string objName;

	string name;

	IloModel model;
	IloCplex cplex;
	IloEnv	 env;

	string dir = "../spInput";
};

void parseCmdLine(int argc, char *argv[], string *inputDir, string *probName, string *outputDir);

int createSGPFInstance(string inputDir, string outputDir);
int createSSNInstance(string inputDir, string outputDir);
int createTransshipInstance(string inputDir, string outputDir);
int createRaoInstance(string inputDir, string outputDir, int argc, char* argv[]);

#endif /* WRITER_HPP_ */
