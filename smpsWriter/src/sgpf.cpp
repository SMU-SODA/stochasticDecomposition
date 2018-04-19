/*
 * sgpf.cpp
 *
 *  Created on: Jan 28, 2018
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

/* The problem first appeared in
 *
 * "SG-Portfolio test problems for stochastic multistage linear programming" by Frauendorfer, K., Hartel, H.,
 * Reiff, M.F. and Schurle, M.
 *
 */

#include "writer.hpp"

class SGPFdata {

public:
	int numMaturities, numAdq, numStd;
	vector<double> ret, initVol;
	vector<int> stdMaturities;
	vector<vector<double>> retBorrow, retLend;

	vector<double> initRet;
};

void defineSGPFdata(SMPSmodel &sgpf, SGPFdata &data);
int createSGPFcor(SMPSmodel &sgpf, SGPFdata &data);
int createSGPFtim(SMPSmodel sgpf);
int createSGPFstoc (SMPSmodel sgpf, SGPFdata data);

int createSGPFInstance(string inputDir, string outputDir) {
	SMPSmodel sgpf;
	SGPFdata data;
	char srcFile[NAMESIZE], destnFile[NAMESIZE];

	/* Setup model parameters */
	defineSGPFdata(sgpf, data);

	/* Generate cor file */
	createSGPFcor(sgpf, data);

	/* Generate tim file */
	createSGPFtim(sgpf);

	/* Generate stoch file */
	createSGPFstoc(sgpf, data);

	/* Move the files to the spInput folder */
	sprintf(srcFile, "sgpf%dy%d.mps", sgpf.numPeriods, sgpf.numStages);
	sprintf(destnFile, "sgpf%dy%d.cor", sgpf.numPeriods, sgpf.numStages);
	rename(srcFile, destnFile);

	sprintf(destnFile, "mkdir %ssgpf%dy%d/", outputDir.c_str(), sgpf.numPeriods, sgpf.numStages);
	system(destnFile);

	sprintf(destnFile, "%ssgpf%dy%d/", outputDir.c_str(), sgpf.numPeriods, sgpf.numStages);
	sprintf(srcFile, "mv sgpf%dy%d.* %s", sgpf.numPeriods, sgpf.numStages, destnFile);
	system(srcFile);

	return 0;
}//END createSGPFInstance()

int createSGPFcor(SMPSmodel &sgpf, SGPFdata &data) {
	char elemName[NAMESIZE];

	try {
		IloEnv   env;
		sprintf(elemName, "sgpf%dy%d", sgpf.numPeriods, sgpf.numStages);
		IloModel model(env, elemName);

		/**************** Decision variables *****************/
		/* assign aliases */
		IloArray<IloNumVarArray> borrow(env, sgpf.numPeriods);
		IloArray<IloNumVarArray> lend(env, sgpf.numPeriods);
		IloArray<IloNumVarArray> volume(env, sgpf.numPeriods);
		IloNumVarArray totalVolume (env, sgpf.numPeriods, 0, IloInfinity);

		/* Load the decision variables onto the solver */
		for ( int t = 0; t < sgpf.numPeriods; t++ ) {
			borrow[t] = IloNumVarArray (env, data.numStd, 0, IloInfinity);
			lend[t] = IloNumVarArray (env, data.numStd, 0, IloInfinity);
			volume[t] = IloNumVarArray (env, data.numMaturities, 0, IloInfinity);

			for ( int i = 0; i < data.numStd; i++ ) {
				/* Volume of maturity borrowed */
				sprintf(elemName, "borrow[%d][%d]", t, data.stdMaturities[i]+1);
				borrow[t][i].setName(elemName); model.add(borrow[t][i]);
				if ( i == 0 )
					sgpf.timCols.push_back(elemName);
				if ( t != 0 ) {
					sgpf.stocCols[t-1].push_back(elemName);
					data.initRet.push_back(data.retBorrow[0][i]);
				}

				/* Volume of maturity lent */
				sprintf(elemName, "lend[%d][%d]", t, data.stdMaturities[i]+1);
				lend[t][i].setName(elemName); model.add(lend[t][i]);
				if ( t != 0 ) {
					sgpf.stocCols[t-1].push_back(elemName);
					data.initRet.push_back(data.retLend[0][i]);
				}
			}

			for ( int i = 0; i < data.numMaturities; i++ ) {
				/* State/volume of the maturity. */
				sprintf(elemName, "volume[%d][%d]", t, i+1);
				volume[t][i].setName(elemName); model.add(volume[t][i]);
			}

			/* State of the portfolio */
			sprintf(elemName, "totalVolume[%d]", t);
			totalVolume[t].setName(elemName); model.add(totalVolume[t]);
			if ( t != 0 ) {
				sgpf.stocCols[t-1].push_back(elemName);
				data.initRet.push_back(data.ret[0]);
			}
		}

		/* Objective function */
		IloObjective obj = IloMinimize(env, 0.0);
		IloExpr netReturn(env);
		for ( int t = 0; t < sgpf.numPeriods; t++) {
			for ( int i = 0; i < data.numStd; i++ ) {
				netReturn += (data.retBorrow[t][i]*borrow[t][i] - data.retLend[t][i]*lend[t][i]);
			}
			netReturn += data.ret[t]*totalVolume[t];
		}
		obj.setExpr(netReturn); model.add(obj); netReturn.end(); sgpf.objName = "obj";

		/* Constraints */
		for ( int t = 0; t < sgpf.numPeriods; t++ ) {
			/* a. State dynamics for a standard maturity */
			for ( int i = 0; i < data.numMaturities; i++ ) {
				IloExpr expr (env);
				sprintf(elemName, "dynamics[%d][%d]", t, i+1);

				vector<int>::iterator it = find(data.stdMaturities.begin(), data.stdMaturities.end(), i);
				if ( it !=  data.stdMaturities.end() ) {
					int j = distance(data.stdMaturities.begin(), it);
					if ( t != 0 && (i+1) != data.numMaturities ) {
						/* The bonds are one year closer to maturity, therefore a _(i+1)_ bond in previous time is now a _i_ time periods
						 * away from maturity. */
						expr = volume[t][i] - volume[t-1][i+1] - borrow[t][j] + lend[t][j];
						IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
					}
					else {
						expr = volume[t][i] - borrow[t][j] + lend[t][j];
						IloConstraint c(expr == data.initVol[i+1]); c.setName(elemName); model.add(c);
					}
				}
				else {
					if ( t != 0) {
						expr = volume[t][i] - volume[t-1][i+1];
						IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
					}
					else {
						expr = volume[t][i];
						IloConstraint c(expr == data.initVol[i+1]); c.setName(elemName); model.add(c);
					}
				}

				if ( i == 0 )
					sgpf.timRows.push_back(elemName);
			}

			/* b. State/volume of the portfolio */
			{
				IloExpr expr (env);
				sprintf(elemName, "pfState[%d]", t);

				expr = totalVolume[t];
				for ( int i = 0; i < data.numMaturities; i++ ) {
					expr -= volume[t][i];
				}
				IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
			}

			/* c. Change in portfolio volume change */
			{
				IloExpr expr (env);
				sprintf(elemName, "volChange[%d]", t);

				double initTotalVol = 0.0;
				expr = totalVolume[t];
				if ( t == 0 ) {
					for ( int i = 0; i < data.numMaturities; i++ ) {
						initTotalVol += data.initVol[i];
					}
				}
				else {
					sgpf.stocRows[t-1].push_back(elemName);
					expr -= totalVolume[t-1];
				}
				IloConstraint c(expr == initTotalVol); c.setName(elemName); model.add(c);
			}

			/* e. Adequacy inequality */
			{
				IloExpr expr (env);
				sprintf(elemName, "adequacy[%d]", t);

				for ( int i = 0; i < data.numMaturities; i++ ) {
					vector<int>::iterator it = find(data.stdMaturities.begin(), data.stdMaturities.end(), i);
					if ( it !=  data.stdMaturities.end() ) {
						int j = distance(data.stdMaturities.begin(), it);
						expr += borrow[t][j];
					}
				}
				double adqTotalVol = 0.0;

				if ( t != 0 ) {
					for ( int i = 0; i < data.numAdq; i++ ) {
						expr -= volume[t-1][i];
					}
					sgpf.stocRows[t-1].push_back(elemName);
				}
				else {
					for ( int i = 0; i < data.numAdq; i++ ) {
						adqTotalVol += data.initVol[i];
					}
				}
				IloConstraint c(expr <= adqTotalVol); c.setName(elemName); model.add(c);
			}
		}

		IloCplex cplex(model);
		sprintf(elemName, "sgpf%dy%d.mps", sgpf.numPeriods, sgpf.numStages);
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
}//END createSGPF()

int createSGPFtim(SMPSmodel sgpf) {
	ofstream tFile;
	char fName[NAMESIZE], line[BLOCKSIZE];

	sprintf(fName, "sgpf%dy%d.tim", sgpf.numPeriods, sgpf.numStages);
	tFile.open(fName);
	tFile << "TIME   sgpf" << sgpf.numPeriods << "y" << sgpf.numStages << endl;

	tFile << "PERIODS" << endl;
	for (int t = 0; t < sgpf.numStages; t++ ) {
		sprintf(line, "    %18s%18s\tStage%2d\n", sgpf.timCols[t].c_str(), sgpf.timRows[t].c_str(), t); tFile << line;
	}
	tFile << "ENDATA" << endl;
	tFile.close();

	return 0;
}//END createSGPFtim()

int createSGPFstoc (SMPSmodel sgpf, SGPFdata data) {
	ofstream sFile;
	char fName[NAMESIZE], line[BLOCKSIZE];
	double vals;

	vals = 0.0;
	for ( int i = 0; i < data.numMaturities; i++ ) {
		vals += data.initVol[i];
	}

	sprintf(fName, "sgpf%dy%d.sto", sgpf.numPeriods, sgpf.numStages); sFile.open(fName);


	sprintf(line, "%-14ssgpf%dy%d\n", "STOCH", sgpf.numPeriods, sgpf.numStages); sFile << line;
	sFile << "INDEP          NORMAL" << endl;

	/* Error terms */
	for (int c = 0; c < (int) sgpf.stocCols[0].size() ; c++ ) {
		sprintf(line, "    U%02d%14s%14.1f%14.1f\n", c+1, "LAGGED", 0.0, 1.0); sFile << line;
	}
	for (int c = 0; c < (int) sgpf.stocRows[0].size() ; c++ ) {
		sprintf(line, "    U%02d%14s%14.1f%14.1f\n", (int) (int) sgpf.stocCols[0].size() + c+1, "LAGGED", 0.0, 0.00005*vals); sFile << line;
	}

	/* Establish the linear relationship between random variables of the stochastic process. */
	sprintf(line, "BLOCKS     LINTR\n"); sFile << line;

	/* Random variables in the model, constant terms based on past values. */
	for ( int t = 1; t < sgpf.numPeriods; t++ ) {
		sprintf(line, " BL PERIOD%d STAGE1\n", t); sFile << line;
		/* Root stage columns and rows do not have uncertain elements, hence index-0 for stocCols/stocRows corresponds
		 * to second stage. */
		for (int c = 0; c < (int) sgpf.stocCols[t-1].size(); c++ ) {
			if ( t == 1 )
				sprintf(line, "    %14s%14s%14f\n", sgpf.stocCols[t-1][c].c_str(), sgpf.objName.c_str(), data.initRet[c]);
			else
				sprintf(line, "    %14s%14s%14f\n", sgpf.stocCols[t-1][c].c_str(), sgpf.objName.c_str(), 0.0);
			sFile << line;
		}
		for ( int i = 0; i < data.numMaturities; i++ ) { /* Computing the initial volume to be used for volume change brownian motion */
			vals += data.initVol[i];
		}
		for (int c = 0; c < (int) sgpf.stocRows[t-1].size(); c++ ) {
			sprintf(line, "    %14s%14s%14f\n", "RHS", sgpf.stocRows[t-1][c].c_str(), 0.001*vals); sFile << line;
		}
	}

	/* Linear transformation matrix */
	for (int c = 0; c < (int) sgpf.stocCols[0].size(); c++ ) {
		sprintf(line, " RV U%02d%11s%14s\n", c+1, "LAGGED", "LAG00"); sFile << line;
		sprintf(line, "    %14s%14s%14f\n", sgpf.stocCols[0][c].c_str(), sgpf.objName.c_str(), 1.0); sFile << line;
	}
	for (int c = 0; c < (int) sgpf.stocRows[0].size(); c++ ) {
		sprintf(line, " RV U%02d%11s%14s\n", (int) sgpf.stocCols[0].size() + c+1, "LAGGED", "LAG00"); sFile << line;
		sprintf(line, "    %14s%14s%14f\n", "RHS", sgpf.stocRows[0][c].c_str(), 1.0); sFile << line;
	}

	for ( int t = 2 ; t < 3; t++ ) {
		for (int c = 0; c < (int) sgpf.stocCols[0].size(); c++ ) {
			sprintf(line, " HV %14s%14s%14s\n", sgpf.stocCols[t-2][c].c_str(), sgpf.objName.c_str(), "LAG01"); sFile << line;
			sprintf(line, "    %14s%14s%14f\n", sgpf.stocCols[t-1][c].c_str(), sgpf.objName.c_str(), 1.0); sFile << line;
		}
		for (int c = 0; c < (int) sgpf.stocRows[0].size(); c++ ) {
			sprintf(line, " HV %14s%14s%14s\n", "RHS", sgpf.stocRows[t-2][c].c_str(), "LAG01"); sFile << line;
			sprintf(line, "    %14s%14s%14f\n", "RHS", sgpf.stocRows[t-1][c].c_str(), 1.0); sFile << line;
		}
	}

	sFile << "ENDATA" << endl;
	sFile.close();

	return 0;
}//END createSGPFtim()

void defineSGPFdata(SMPSmodel &sgpf, SGPFdata &data) {
	int parse;

	/* Number of periods to be considered and the number of stages. */
	cout << "Enter model parameters (1)/use default values (0) : ";
	cin >> parse;

	if ( parse ) {
		cout << "Enter the number of periods in the model : ";
		cin >> sgpf.numPeriods;
		cout << "Enter the number of stages               : ";
		cin >> sgpf.numStages;
	}
	else {
		sgpf.numPeriods = 3;
		sgpf.numStages = 2;
	}

	sgpf.stocCols = vector<vector<string>> (sgpf.numPeriods, vector<string> ());
	sgpf.stocRows = vector<vector<string>> (sgpf.numPeriods, vector<string> ());

	data.numMaturities = 60;	/* Total number of maturities */
	data.numAdq = 3;			/* Number of maturities used to ensure adequacy (M in the paper) */
	data.stdMaturities = {1,2,3,6,12,24,36,48,60}; /* list of standard maturities */
	data.numStd = (int) data.stdMaturities.size();

	/* Change list of standard maturities to suit indexing */
	for ( int n = 0; n < data.numStd; n++ ) {
		data.stdMaturities[n] -= 1;
	}

	/* Return rate for return, borrowing, lending */
	data.ret 		= vector<double> (sgpf.numPeriods);
	data.retBorrow 	= vector< vector<double>> (sgpf.numPeriods, vector<double> (data.numStd) );
	data.retLend 	= vector< vector<double>> (sgpf.numPeriods, vector<double> (data.numStd) );
	data.initVol 	= vector<double> (data.numMaturities);

	default_random_engine generator{static_cast<long unsigned int>(3554548844580680)};
	normal_distribution <double> distribution(0.0,1.0);

	/* Take care of the boundary conditions. Remove the first 200 (ad-hoc choice) observations. */
	double dataInit = 0;
	{
		int t = 0;
		while ( t < 200 ) {
			dataInit += distribution(generator);
			t++;
		}
	}

	/* Returns */
	for ( int t = 0; t < sgpf.numPeriods; t++ ) {
		if ( t == 0 ) {
			data.ret[t] = dataInit + distribution(generator);
 		}
		else {
			data.ret[t] = data.ret[t-1] + distribution(generator);
		}
	}

	/* Borrowing rate */
	for ( int n = 0; n < data.numStd; n++ ) {
		for ( int t = 0; t < sgpf.numPeriods; t++ ) {
			double number;
			do {
				number = distribution(generator);
			}
			while ( t == 0 && ((data.retBorrow[t][n] = number) < TOLERANCE) );

			if ( t != 0 )
				while ( (data.retBorrow[t][n] = data.retBorrow[t-1][n] + number) < TOLERANCE )
					number = distribution(generator);
		}
	}

	/* Lending rate */
	for ( int n = 0; n < data.numStd; n++ ) {
		for ( int t = 0; t < sgpf.numPeriods; t++ ) {
			double number;
			do {
				number = distribution(generator);
			}
			while ( t == 0 && ((data.retLend[t][n] = number) < TOLERANCE) );

			if ( t != 0 )
				while ( (data.retLend[t][n] = data.retBorrow[t-1][n] + number) < TOLERANCE )
					number = distribution(generator);
		}
	}

	/* Initial volume */
	data.initVol[0] = 200000; data.initVol[1] = 120000; data.initVol[2] = 80000;

	return;
}//END defineData()
