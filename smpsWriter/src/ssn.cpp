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

	vector<vector<double>> vals, probs;					// Values and probabilities of stochastic elements in the problem (demand).
	vector<double> meanDemand;

	int numGroups;										// Clustering all the source-destination pairs which have these many routes */
	vector<vector<int>> groupMembers;

	SSNdata(SMPSmodel &ssn, string &inputDir);
	~SSNdata();
};

int createSSNcor(SMPSmodel &ssn, SSNdata &data);
int createSSNtim(SMPSmodel ssn, SSNdata data);
int createSSNstoc (SMPSmodel ssn, SSNdata data);

int createSSNInstance(string inputDir, string outputDir) {
	SMPSmodel ssn;
	char srcFile[NAMESIZE], destnFile[NAMESIZE];
	int parse;

	/* Read the data relevant to the problem */
	SSNdata data(ssn, inputDir);

	/* Number of periods to be considered and the number of stages. */
	cout << "Enter model parameters (1)/use default values (0) : ";
	cin >> parse;

	if ( parse == 1 ) {
		cout << "Enter the number of groups : ";
		cin >> data.numGroups;
	}
	data.groupMembers = vector<vector<int>> (data.numGroups);
	for ( int i = 0; i < (int) data.numSDpairs; i++ ) {
		int idx = data.meanDemand[i]/10 < (data.numGroups-1)?
				data.meanDemand[i]/10:(data.numGroups-1);
		data.groupMembers[idx].push_back(i);
	}

	/* Create the SMPS file */
	createSSNcor(ssn, data);
	createSSNtim(ssn, data);
	createSSNstoc(ssn, data);

	/* Move the files to the spInput folder */
	sprintf(destnFile, "ssn_rc%d.cor", data.numGroups);
	sprintf(srcFile, "ssn_rc%d.mps", data.numGroups);
	rename(srcFile, destnFile);

	sprintf(destnFile, "mkdir %sssn_rc%d/", outputDir.c_str(), data.numGroups);
	system(destnFile);

	sprintf(destnFile, "%sssn_rc%d/", outputDir.c_str(), data.numGroups);
	sprintf(srcFile, "mv ssn_rc%d.* %s", data.numGroups, destnFile);
	system(srcFile);

	printf("Successfully wrote the SMPS files for %s%d to the output directory '%s'.\n", "ssn_rc", data.numGroups, outputDir.c_str());

	return 0;
}//END createSSNInstance()

int createSSNcor(SMPSmodel &ssn, SSNdata &data) {
	char elemName[NAMESIZE];

	try {
		IloEnv env;
		sprintf(elemName, "ssn_rc%d", data.numGroups);
		IloModel model(env, elemName);

		/**************** Decision variables *****************/
		IloNumVarArray capAdd(env, data.numLinks, 0, IloInfinity, ILOFLOAT);
		IloNumVarArray dropCall(env, data.numSDpairs, 0, IloInfinity, ILOFLOAT);
		IloArray<IloNumVarArray> flow(env, data.numSDpairs);
		IloNumVarArray group(env, data.numGroups, 0, IloInfinity, ILOFLOAT);

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

		/* Grouping of source-destination pairs */
		for (int i = 0; i < data.numGroups; i++ ) {
			sprintf(elemName, "group[%d]", i+1);
			group[i].setName(elemName); model.add(group[i]);
			ssn.stocCols[0].push_back(elemName);
		}

		/******************* Constraints *********************/
		/* Budget limits on expansion */
		IloExpr expr (env);
		sprintf(elemName, "budget");
		for ( int i = 0; i < data.numLinks; i++ ) {
			expr += capAdd[i];
		}
		IloConstraint c(expr <= data.budget); c.setName(elemName); model.add(c);
		ssn.timRows.push_back(elemName);

		/* Flow balance for each source-destination pair */
		for ( int i = 0; i < data.numSDpairs; i++ ) {
			IloExpr expr(env);
			sprintf(elemName, "demand[%d]", i);

			for ( int r = 0; r < (int) data.routes4SDpair[i].size(); r++ ) {
				expr += flow[i][r];
			}
			expr += dropCall[i];
			IloConstraint c(expr == data.demand[i]); c.setName(elemName); model.add(c);

			ssn.stocRows[0].push_back(elemName);
			if ( i == 0 )
				ssn.timRows.push_back(elemName);
		}

		/* Capacity limits (existing and added) */
		for ( int j = 0; j < data.numLinks; j++ ) {
			IloExpr expr(env);
			sprintf(elemName, "linkCap[%d]", j);

			for ( int r = 0; r < (int) data.routesOnLinks[j].size(); r++ ) {
				int i = 0; vector<int>::iterator it;
				while ((it = find(data.routes4SDpair[i].begin(), data.routes4SDpair[i].end(), data.routesOnLinks[j][r]))
						== data.routes4SDpair[i].end()) {
					i++;
					if ( i == (int) data.routes4SDpair.size() )
						perror("Could not match the route to the any source-destination pair.\n");
				}
				ptrdiff_t pos = it - data.routes4SDpair[i].begin();
				expr += flow[i][pos];
			}
			expr -= capAdd[j];

			IloConstraint c(expr <= data.capacities[j]); c.setName(elemName); model.add(c);
		}

		/* Group the source-destination pairs based on given criterion */
		for ( int g = 0; g < data.numGroups; g++ ) {
			IloExpr expr(env);
			sprintf(elemName, "set[%d]", g);
			expr = group[g];
			for ( int i = 0; i < (int) data.groupMembers[g].size(); i++ ) {
				expr -= dropCall[data.groupMembers[g][i]];
			}
			IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
		}

		/************** Objective Function *******************/
		IloObjective obj = IloMinimize(env, 0.0);
		IloExpr losses(env);
		for ( int i = 0; i < data.numGroups; i++ ) {
			losses += group[i];
		}
		obj.setExpr(losses); model.add(obj); losses.end(); ssn.objName = "obj";

		IloCplex cplex(model);
		sprintf(elemName, "ssn_rc%d.mps", data.numGroups);
		cplex.exportModel(elemName);
		sprintf(elemName, "ssn_rc%d.lp", data.numGroups);
		cplex.exportModel(elemName);
		env.end();
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}

	return 0;
}//END createSSNcor()

int createSSNtim(SMPSmodel ssn, SSNdata data) {
	ofstream tFile;
	char fName[NAMESIZE], line[BLOCKSIZE];

	sprintf(fName, "ssn_rc%d.tim", data.numGroups);
	tFile.open(fName);
	tFile << "TIME   ssn_rc" << endl;

	tFile << "PERIODS" << endl;
	for (int t = 0; t < ssn.numStages; t++ ) {
		sprintf(line, "    %18s%18s\tStage%02d\n", ssn.timCols[t].c_str(), ssn.timRows[t].c_str(), t); tFile << line;
	}
	tFile << "ENDATA" << endl;
	tFile.close();

	return 0;
}//END createSSNtim()

int createSSNstoc (SMPSmodel ssn, SSNdata data) {
	ofstream sFile;
	char fName[NAMESIZE], line[BLOCKSIZE];

	sprintf(fName, "ssn_rc%d.sto", data.numGroups); sFile.open(fName);

	sprintf(line, "%-14sssn_rc\n", "STOCH"); sFile << line;
	sFile << "INDEP          DISCRETE" << endl;

	/* Error terms */
	for (int i = 0; i < (int) ssn.stocRows[0].size() ; i++ ) {
		for ( int j = 0; j < (int) data.vals[i].size(); j++ ) {
			sprintf(line, "    RHS%14s%14.5f%14.5f\n", ssn.stocRows[0][i].c_str(), data.vals[i][j], data.probs[i][j]); sFile << line;
		}
	}
	for (int i = 0; i < data.numGroups; i++ ) {
		double vals = 10/(i+1);
		sprintf(line, "    %s%14s%14.5f%14.5f\n", ssn.stocCols[0][i].c_str(), ssn.objName.c_str(), 0.2*vals, 0.1); sFile << line;
		sprintf(line, "    %s%14s%14.5f%14.5f\n", ssn.stocCols[0][i].c_str(), ssn.objName.c_str(), 0.5*vals, 0.2); sFile << line;
		sprintf(line, "    %s%14s%14.5f%14.5f\n", ssn.stocCols[0][i].c_str(), ssn.objName.c_str(), 1.0*vals, 0.4); sFile << line;
		sprintf(line, "    %s%14s%14.5f%14.5f\n", ssn.stocCols[0][i].c_str(), ssn.objName.c_str(), 1.5*vals, 0.2); sFile << line;
		sprintf(line, "    %s%14s%14.5f%14.5f\n", ssn.stocCols[0][i].c_str(), ssn.objName.c_str(), 2.0*vals, 0.1); sFile << line;
	}

	sFile << "ENDATA" << endl;
	sFile.close();

	return 0;
}//END createSGPFtim()

SSNdata::SSNdata(SMPSmodel &ssn, string &inputDir) {
	vector<string> rvRow;

	ssn.numPeriods = 1;
	ssn.numStages = 2;

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
						currentRoutes.push_back(atoi(tokens[i].c_str())-1);
					}
					routesOnLinks.push_back(currentRoutes);
					j++;
				}
				else if ( section == 2 )  {
					demand.push_back(atof(tokens[0].c_str()));
					for ( int i = 1; i < (int) tokens.size(); i++ ) {
						currentRoutes.push_back(atoi(tokens[i].c_str())-1);
					}
					routes4SDpair.push_back(currentRoutes);
				}
				else if ( section == 3 ) {
					vector<string>::iterator it; ptrdiff_t pos;
					if ( (it = find(rvRow.begin(), rvRow.end(), tokens[1])) == rvRow.end() ) {
						pos = (int) rvRow.size();
						rvRow.push_back(tokens[1]);
						if ( pos == 0 ) {
							vals = vector<vector<double>> (numSDpairs);
							probs = vector<vector<double>> (numSDpairs);
						}
						vals[pos] = vector<double> ();
						probs[pos] = vector<double> ();
					}
					else {
						pos = it - rvRow.begin();
					}
					vals[pos].push_back(atof(tokens[2].c_str()));
					probs[pos].push_back(atof(tokens[3].c_str()));
				}
			}
			else {
				section++;
				j = 0;
			}
		}
	}
	else
		perror("Failed to open the data file for SSN.\n");
	ssnDataFile.close();

	meanDemand = vector<double> ((int) vals.size());
	/* Compute the mean demand */
	for ( int i = 0; i < (int) vals.size(); i++ ) {
		meanDemand[i] = 0;
		for ( int j = 0; j < (int) vals[i].size(); j++ ) {
			meanDemand[i] += vals[i][j]*probs[i][j];
		}
	}

	/* prepare to store the columns and rows with random elements in them. */
	ssn.stocCols = vector<vector<string>> (ssn.numPeriods, vector<string> ());
	ssn.stocRows = vector<vector<string>> (ssn.numPeriods, vector<string> ());

	numGroups = 1;
}

SSNdata::~SSNdata() {

}
