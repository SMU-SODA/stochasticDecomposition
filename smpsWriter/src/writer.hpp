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
#include <algorithm>

#include <ilcplex/ilocplex.h>

using namespace std;

#define NAMESIZE 32

class SMPSmodel {

public:
	SMPSmodel();
	~SMPSmodel();

	int numStages, numPeriods;
	vector<string> timCols, timRows, stocRows, stocCols;
	string objName;
};

void parseCmdLine(int argc, char *argv[], string *probName);

int createSGPFInstance();
int createSGPFcor(SMPSmodel &sgpf);
int createSGPFtim(SMPSmodel sgpf);
void defineSGPFData();

#endif /* WRITER_HPP_ */
