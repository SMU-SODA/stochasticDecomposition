/*
 * ssn.cpp
 *
 *  Created on: Apr 11, 2018
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "writer.hpp"

class SSNdata{

public:
	int numLinks, numSDpairs;							// Number of links and source-destination pairs
	double budget;										// Budget for capacity expansion
	vector<vector<int>> routesOnLinks, routes4SDpair;	// All routes on a link and all routes between source-destination pair
	vector<double> capacities, demand;					// Existing capacities and demand for each source-destination pair

	SSNdata(string &inputDir);
	~SSNdata();
};

int createSSNcor(SMPSmodel &ssn, SSNdata &data);

int createSSNInstance(string inputDir, string outputDir) {
	SMPSmodel ssn;
	char srcFile[NAMESIZE], destnFile[NAMESIZE];

	SSNdata data(inputDir);

	/* Create the core file */
	createSSNcor(ssn, data);

	/* Move the files to the spInput folder */
	rename("ssn.mps", "ssn.cor");

	sprintf(destnFile, "mkdir %sssn/", outputDir.c_str());
	system(destnFile);

	sprintf(destnFile, "%sssn/", outputDir.c_str());
	sprintf(srcFile, "mv ssn.* %s", destnFile);
	system(srcFile);

	return 0;
}//END createSSNInstance()

int createSSNcor(SMPSmodel &ssn, SSNdata &data) {
	char elemName[NAMESIZE];

	try {
		IloEnv env;
		sprintf(elemName, "ssn");
		IloModel model(env, elemName);

		/**************** Decision variables *****************/
		IloNumVarArray capAdd(env, data.numLinks);
		IloNumVarArray dropCall(env, data.numSDpairs);
		IloArray<IloNumVarArray> flow(env, data.numSDpairs);

		/* Capacity extension */
		for ( int j = 0; j < data.numLinks; j++ ) {
			sprintf(elemName, "capAdd[%d]", j);
			capAdd[j].setName(elemName); model.add(capAdd[j]);
			if ( j == 0 ) {
				ssn.timCols.push_back(elemName);
			}
		}

		/* Flow on a route */
		for ( int i = 0; i < data.numSDpairs; i++ ) {
			flow[i] = IloNumVarArray(env, (int) data.routes4SDpair[i].size(), 0, IloInfinity, ILOFLOAT);
			for ( int r = 0; r < (int) data.routes4SDpair[i].size(); r++ ) {
				sprintf(elemName, "flow[%d][%d]", i, r);
				flow[i][r].setName(elemName); model.add(flow[i][r]);
				if ( i == 0 && r == 0 ) {
					ssn.timCols.push_back(elemName);
				}
			}
		}

		/* Dropped call/demand */
		for ( int i = 0; i < data.numSDpairs; i++ ) {
			sprintf(elemName, "dropCall[%d]", i);
			dropCall[i].setName(elemName); model.add(dropCall[i]);
		}

		/******************* Constraints *********************/
		/* Budget limits on expansion */
		IloExpr expr (env);
		sprintf(elemName, "budget");
		for ( int i = 0; i < data.numLinks; i++ ) {
			expr += capAdd[i];
			IloConstraint c(expr <= data.budget); c.setName(elemName); model.add(c);
		}

		/* Capacity limits (existing and added) */
		for ( int j = 0; j < data.numLinks; j++ ) {
			IloExpr expr(env);
			sprintf(elemName, "linkCap[%d]", j);

			for ( int i = 0; i < data.numSDpairs; i++ ) {
				for ( int r = 0; r < (int) data.routes4SDpair.size(); r++ ) {
					expr += flow[i][r];
				}
			}
			expr -= capAdd[j];

			IloConstraint c(expr <= data.capacities[j]); c.setName(elemName); model.add(c);
		}

		/* Flow balance for each source-destination pair */
		for ( int i = 0; i < data.numSDpairs; i++ ) {
			IloExpr expr(env);
			sprintf(elemName, "demand[%d]", i);

			for ( int r = 0; r < (int) data.routes4SDpair.size(); r++ ) {
				expr += flow[i][r];
			}
			expr += dropCall[i];
			IloConstraint c(expr == data.demand[i]); c.setName(elemName); model.add(c);
		}

		/************** Objective Function *******************/
		IloObjective obj = IloMinimize(env, 0.0);
		IloExpr losses(env);
		for ( int i = 0; i < data.numSDpairs; i++ ) {
			losses += dropCall[i];
		}
		obj.setExpr(losses); model.add(obj); losses.end(); ssn.objName = "obj";

	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}

	return 0;
}//END createSSNcor()

SSNdata::SSNdata(string &inputDir) {

	numLinks = 89;
	numSDpairs = 86;

	budget = 1008;

	/* Read the data file for SSN to extract problem information */
	ifstream ssnDataFile;
	ssnDataFile.open((inputDir + "ssnData.txt").c_str());
	if (ssnDataFile.is_open()) {
		string line;
		int section = 0, j;
		while ( getline(ssnDataFile, line) ){
			if ( line[0] != '#' ) {
				istringstream iss(line);
				vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};
				vector<int> currentRoutes;

				if ( section == 1 ) {
					capacities.push_back(atoi(tokens[0].c_str()));
					for ( int i = 1; i < (int) tokens.size(); i++ ) {
						currentRoutes.push_back(atoi(tokens[i].c_str()));
					}
					routesOnLinks.push_back(currentRoutes);
					j++;
				}
				else if ( section == 2 )  {
					demand.push_back(atof(tokens[0].c_str()));
					for ( int i = 1; i < (int) tokens.size(); i++ ) {
						currentRoutes.push_back(atoi(tokens[i].c_str()));
					}
					routes4SDpair.push_back(currentRoutes);
				}
			}
			else {
				section++;
				j = 0;
			}
		}
	}
	ssnDataFile.close();
}

SSNdata::~SSNdata() {

}
