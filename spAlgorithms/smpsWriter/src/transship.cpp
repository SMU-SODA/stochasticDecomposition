/*
 * transship.cpp
 *
 *  Created on: May 18, 2018
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */
#include "writer.hpp"

class TransshipData {

public:
	int numRetailers;

	vector<double> holdingCost, backlogCost, replenishCost;
	vector<vector<double>> transshipCost;

	vector<double> meanDem, sdDem;
};

int createTransshipCor(SMPSmodel &transship, TransshipData &data);
int createTransshipTim(SMPSmodel trans);
int createTransshipStoc (SMPSmodel trans, TransshipData &data);
void defineTransshipData(SMPSmodel &transship, TransshipData &data);

int createTransshipInstance(string inputDir, string outputDir) {
	SMPSmodel trans;
	TransshipData data;
	char srcFile[NAMESIZE], destnFile[NAMESIZE];

	/* Setup model parameters */
	defineTransshipData(trans, data);

	/* Generate cor file */
	createTransshipCor(trans, data);
	createTransshipTim(trans);
	createTransshipStoc(trans, data);

	/* Move the files to the spInput folder */
	sprintf(destnFile, "transship.cor");
	sprintf(srcFile, "transship.mps");
	rename(srcFile, destnFile);

	sprintf(destnFile, "mkdir %stransship/", outputDir.c_str());
	system(destnFile);

	sprintf(destnFile, "%stransship/", outputDir.c_str());
#ifdef _WIN64
	sprintf(srcFile, "move transship.* %s", destnFile);
#else
	sprintf(srcFile, "mv transship.* %s", destnFile);
#endif

	system(srcFile);

	printf("Successfully wrote the SMPS files for %s to the output directory '%s'.\n", "transshipment",
			outputDir.c_str());

	return 0;
}//END createTransshipInstance()

int createTransshipCor(SMPSmodel &trans, TransshipData &data) {
	char elemName[NAMESIZE];

	try {
		IloEnv   env;
		sprintf(elemName, "transship");
		IloModel model(env, elemName);

		/**************** Decision variables *****************/
		trans.timeCols.push_back("orderUp(0)");
		IloNumVarArray orderUp(env, data.numRetailers, 0, IloInfinity);
		orderUp.setNames("orderUp"); model.add(orderUp);

		trans.timeCols.push_back("begEnd(0)");
		IloNumVarArray begEnd(env, data.numRetailers, 0, IloInfinity);
		begEnd.setNames("begEnd"); model.add(begEnd);

		IloNumVarArray repDem(env, data.numRetailers, 0, IloInfinity);
		repDem.setNames("repDem"); model.add(repDem);
		for ( int n = 0; n < data.numRetailers; n++ ) {
			sprintf(elemName, "repDem(%d)", n);
			trans.stocCols[0].push_back(elemName);
			trans.stocRows[0].push_back("obj");
		}

		IloNumVarArray repEnd(env, data.numRetailers, 0, IloInfinity);
		repEnd.setNames("repEnd"); model.add(repEnd);

		IloNumVarArray demVar(env, data.numRetailers, 0, IloInfinity);
		demVar.setNames("demVar"); model.add(demVar);

		IloArray<IloNumVarArray> begDem(env, data.numRetailers);

		/* Load the decision variables onto the solver */
		for ( int n = 0; n < data.numRetailers; n++ ) {
			begDem[n] = IloNumVarArray(env, data.numRetailers, 0, IloInfinity);
			sprintf(elemName, "begDem(%d)", n);
			begDem[n].setNames(elemName);
		}

		/**************** Objective function *****************/
		IloObjective obj = IloMinimize(env, 0.0);
		IloExpr totalCost(env);
		for ( int n = 0; n < data.numRetailers; n++ ) {
			totalCost += data.holdingCost[n]*begEnd[n] + data.backlogCost[n]*repDem[n];
			for ( int m = 0; m < data.numRetailers; m++ ) {
				totalCost += (data.transshipCost[n][m] + data.replenishCost[n] - data.replenishCost[m])*begDem[n][m];
			}
		}
		obj.setExpr(totalCost); model.add(obj); totalCost.end(); trans.objName = "obj";
		trans.timeRows.push_back("obj");

		/**************** Constraints *****************/
		/* Initial inventory at retailer */
		for ( int n = 0; n < data.numRetailers; n++ ) {
			IloExpr expr (env);
			sprintf(elemName, "initInv(%d)", n);
			if ( n == 0) {
				trans.timeRows.push_back(elemName);
			}

			expr = - orderUp[n] + begDem[n][n] + begEnd[n];
			for ( int m = 0; m < data.numRetailers; m++ ) {
				if ( m != n ) {
					expr += begDem[n][m];
				}
			}
			IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
		}

		/* Sink node at the retailer */
		for ( int n = 0; n < data.numRetailers; n++ ) {
			IloExpr expr (env);
			sprintf(elemName, "retailDem(%d)", n);

			expr = begDem[n][n] + repDem[n] - demVar[n];
			for ( int m = 0; m < data.numRetailers; m++ ) {
				if ( m != n ) {
					expr += begDem[m][n];
				}
			}
			IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
		}

		/* Terminal inventory at the retailer */
		for ( int n = 0; n < data.numRetailers; n++ ) {
			IloExpr expr (env);
			sprintf(elemName, "finalInv(%d)", n);

			expr = begEnd[n] + repEnd[n] - orderUp[n];
			IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
		}

		/* Replenishment nodes */
		for ( int n = 0; n < data.numRetailers; n++ ) {
			IloExpr expr (env);
			sprintf(elemName, "repVal(%d)", n);

			expr += (repDem[n] + repEnd[n] - demVar[n]);
			IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
		}

		/* Auxiliary constraint to conform with SMPS format */
		for ( int n = 0; n < data.numRetailers; n++ ) {
			IloExpr expr(env);
			sprintf(elemName, "dummy(%d)", n);
			trans.stocRows[0].push_back(elemName);
			trans.stocCols[0].push_back("RHS");

			expr = demVar[n];
			IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
		}

		IloCplex cplex(model);
		sprintf(elemName, "transship.mps");
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
}//END createTransshipCor()

int createTransshipTim(SMPSmodel trans) {
	ofstream tFile;
	char fName[NAMESIZE], line[BLOCKSIZE];

	sprintf(fName, "transship.tim");
	tFile.open(fName);
	tFile << "TIME   transship" << endl;

	tFile << "PERIODS" << endl;
	for (int t = 0; t < trans.numStages; t++ ) {
		sprintf(line, "    %18s%18s\tStage%2d\n", trans.timeCols[t].c_str(), trans.timeRows[t].c_str(), t); tFile << line;
	}
	tFile << "ENDATA" << endl;
	tFile.close();

	return 0;
}//END createSGPFtim()

int createTransshipStoc (SMPSmodel trans, TransshipData &data) {
	ofstream sFile;
	char fName[NAMESIZE], line[BLOCKSIZE];

	sprintf(fName, "transship.sto"); sFile.open(fName);

	sprintf(line, "%-14s%-14s\n", "STOCH", "transship"); sFile << line;
	sFile << "INDEP          NORMAL" << endl;

	int idx = 0;
	for (int c = 0; c < (int) trans.stocCols[0].size() ; c++ ) {
		if ( trans.stocCols[0][c] == "RHS" ) {
			sprintf(line, "    %s%14s%14.5f%14.5f\n", trans.stocCols[0][c].c_str(),
					trans.stocRows[0][c].c_str(), data.meanDem[idx], data.sdDem[idx]);
			sFile << line;
			idx++;
		}
	}

	idx = 0;
	for (int c = 0; c < (int) trans.stocCols[0].size() ; c++ ) {
		if (trans.stocCols[0][c] != "RHS" ) {
			sprintf(line, "    %s%14s%14.5f%14.5f\n", trans.stocCols[0][c].c_str(),
					trans.stocRows[0][c].c_str(), data.backlogCost[idx], 0.2*data.backlogCost[idx]);
			idx++;
			sFile << line;
		}
	}

	sFile << "ENDATA" << endl;
	sFile.close();

	return 0;
}

//for ( int i = 0; i < (int) trans.stocRows[0].size(); i++ ) {
//}

void defineTransshipData(SMPSmodel &trans, TransshipData &data) {

	data.numRetailers = 7;

	data.holdingCost = {1.00, 1.05, 1.1, 1.15, 1.20, 1.25, 1.3};
	data.backlogCost = {4.00, 4.20, 4.40, 4.60, 4.80, 5.00, 5.20};

	data.transshipCost = vector<vector<double>> (data.numRetailers, vector<double> (data.numRetailers, 2.00));
	data.replenishCost = vector<double> (data.numRetailers, 0.5);

	data.meanDem = {100,200,150,170,180,170,170};
	data.sdDem = {20,50,30,50,40,30,50};

	/* prepare to store the columns and rows with random elements in them. */
	trans.stocCols = vector<vector<string>> (1, vector<string> ());
	trans.stocRows = vector<vector<string>> (1, vector<string> ());


}//END defineTransshipData()
