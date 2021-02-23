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
	vector<vector<double>> retBorrow;

	double delta;

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
		IloArray<IloNumVarArray> volChg(env, sgpf.numPeriods);
		IloArray<IloNumVarArray> lend(env, sgpf.numPeriods);
		IloArray<IloNumVarArray> vol(env, sgpf.numPeriods);
		IloNumVarArray totalVol (env, sgpf.numPeriods, -IloInfinity, IloInfinity);

		/* Load the decision variables onto the solver */
		for ( int t = 0; t < sgpf.numPeriods; t++ ) {
			volChg[t] = IloNumVarArray (env, data.numStd, -IloInfinity, IloInfinity);
			lend[t] = IloNumVarArray (env, data.numStd, 0, IloInfinity);
			vol[t] = IloNumVarArray (env, data.numMaturities, 0, IloInfinity);

			for ( int i = 0; i < data.numStd; i++ ) {
				/* Volume of maturity borrowed */
				sprintf(elemName, "volChg[%d][%d]", t, data.stdMaturities[i]+1);
				volChg[t][i].setName(elemName); model.add(volChg[t][i]);
				if ( i == 0 )
					sgpf.timeCols.push_back(elemName);
				if ( t != 0 ) {
					sgpf.stocCols[t-1].push_back(elemName);
					sgpf.stocRows[t-1].push_back("obj");
					data.initRet.push_back(data.retBorrow[0][i]);
				}

				/* Volume of maturity lent */
				sprintf(elemName, "lend[%d][%d]", t, data.stdMaturities[i]+1);
				lend[t][i].setName(elemName); model.add(lend[t][i]);
			}

			for ( int i = 0; i < data.numMaturities; i++ ) {
				/* State/volume of the maturity. */
				sprintf(elemName, "vol[%d][%d]", t, i+1);
				vol[t][i].setName(elemName); model.add(vol[t][i]);
			}

			/* State of the portfolio */
			sprintf(elemName, "totalVol[%d]", t);
			totalVol[t].setName(elemName); model.add(totalVol[t]);
			if ( t != 0 ) {
				sgpf.stocCols[t-1].push_back(elemName);
				sgpf.stocRows[t-1].push_back("obj");
				data.initRet.push_back(data.ret[0]);
			}
		}

		/* Objective function */
		IloObjective obj = IloMinimize(env, 0.0);
		IloExpr netReturn(env);
		for ( int t = 0; t < sgpf.numPeriods; t++) {
			for ( int i = 0; i < data.numStd; i++ ) {
				netReturn += (data.retBorrow[t][i]*volChg[t][i] + data.delta*lend[t][i]);
			}
			netReturn -= data.ret[t]*totalVol[t];
		}
		obj.setExpr(netReturn); model.add(obj); netReturn.end(); sgpf.objName = "obj";

		/* Constraints */
		for ( int t = 0; t < sgpf.numPeriods; t++ ) {
			/* Evolution of the maturity. */
			for ( int i = 0; i < data.numMaturities; i++ ) {
				IloExpr expr (env);
				sprintf(elemName, "matDyn[%d][%d]", t, i+1);

				vector<int>::iterator it = find(data.stdMaturities.begin(), data.stdMaturities.end(), i);
				if ( it !=  data.stdMaturities.end() ) {
					/* 1a. State dynamics for a standard maturity */
					int j = distance(data.stdMaturities.begin(), it);
					if ( t != 0 ) {
						/* The bonds are one year closer to maturity, therefore a _(i+1)_ bond in previous time is now a _i_ time periods
						 * away from maturity. */
						if ( (i+1) != data.numMaturities ) {
							expr = vol[t][i] - vol[t-1][i+1] - volChg[t][j];
							IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
						}
						else {
							/* For the final maturity (60 months), you can only borrow. There is no money in this bond to lend. */
							expr = vol[t][i] - volChg[t][j];
							IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
						}
					}
					else {
						expr = vol[t][i] - volChg[t][j];
						IloConstraint c(expr == data.initVol[i+1]); c.setName(elemName); model.add(c);
					}
				}
				else {
					/* 1b. State dynamics for a non-standard maturity */
					if ( t != 0) {
						expr = vol[t][i] - vol[t-1][i+1];
						IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
					}
					else {
						expr = vol[t][i];
						IloConstraint c(expr == data.initVol[i+1]); c.setName(elemName); model.add(c);
					}
				}

				if ( i == 0 )
					sgpf.timeRows.push_back(elemName);
			}

			/* 2. State/volume of the portfolio */
			{
				IloExpr expr (env);
				sprintf(elemName, "pfState[%d]", t);

				expr = totalVol[t];
				for ( int i = 0; i < data.numMaturities; i++ ) {
					expr -= vol[t][i];
				}
				IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
			}

			/* 3. Change in portfolio volume change */
			{
				IloExpr expr (env);
				sprintf(elemName, "pfChange[%d]", t);

				double initTotalVol = 0.0;
				if ( t == 0 ) {
					expr = totalVol[t];
					for ( int i = 0; i < data.numMaturities; i++ ) {
						initTotalVol += data.initVol[i];
					}
					/* For the first time period, we allow a maximum of 20% change in portfolio volume. */
					initTotalVol *= 0.8;
				}
				else {
					expr = totalVol[t] - totalVol[t-1];
				}
				IloConstraint c(expr == initTotalVol); c.setName(elemName); model.add(c);

				sgpf.stocRows[t-1].push_back(elemName);
				sgpf.stocCols[t-1].push_back("RHS");
			}

			/* 4. Adequacy inequality */
			{
				IloExpr expr (env);
				sprintf(elemName, "pfAdq[%d]", t);

				for ( int i = 0; i < data.numMaturities; i++ ) {
					vector<int>::iterator it = find(data.stdMaturities.begin(), data.stdMaturities.end(), i);
					if ( it !=  data.stdMaturities.end() ) {
						int j = distance(data.stdMaturities.begin(), it);
						expr += (volChg[t][j] + lend[t][j]);
					}
				}

				double adqTotalVol = 0.0;
				if ( t != 0 ) {
					for ( int i = 0; i < data.numAdq; i++ ) {
						expr -= vol[t-1][i];
					}
					expr -= (totalVol[t] - totalVol[t-1]);
				}
				else {
					expr -= totalVol[t];

					for ( int i = data.numAdq+1; i < data.numMaturities; i++ ) {
						adqTotalVol -= data.initVol[i];
					}
				}
				IloConstraint c(expr <= adqTotalVol); c.setName(elemName); model.add(c);
			}

			/* 5. Relation between borrowing and lending amounts. */
			{
				for ( int i = 0; i < data.numMaturities; i++ ) {
					vector<int>::iterator it = find(data.stdMaturities.begin(), data.stdMaturities.end(), i);
					if ( it !=  data.stdMaturities.end() ) {
						int j = distance(data.stdMaturities.begin(), it);

						IloExpr expr (env);
						sprintf(elemName, "blReln[%d][%d]", t, i+1);

						expr = volChg[t][j] + lend[t][j];

						IloConstraint c(expr >= 0); c.setName(elemName); model.add(c);
					}
				}
			}
		}

		IloCplex cplex(model);
#if 0
		sprintf(elemName, "sgpf%dy%d.mps", sgpf.numPeriods, sgpf.numStages);
		cplex.exportModel(elemName);
#endif

		cplex.solve();
		cout << "solution status = " << cplex.getStatus() << endl;

		cout << "Optimal objective function value = " << cplex.getObjValue() << endl;
		cout << "Optimal solution = " << cplex.getObjValue() << endl;

		cout << "  Total volume = " << cplex.getValue(totalVol[0]) << endl;
		for ( int t = 0; t < sgpf.numPeriods; t++ ) {
			for (int i = 0; i < data.numStd; i++) {
				fprintf(stdout, "Buy[%d][%d] = %lf;\t Sell[%d][%d] = %lf;\t Volume[%d][%d] = %lf\n",
						t, data.stdMaturities[i]+1, cplex.getValue(volChg[t][i] + lend[t][i]),
						t, data.stdMaturities[i]+1, cplex.getValue(lend[t][i]),
						t, data.stdMaturities[i]+1, cplex.getValue(vol[t][data.stdMaturities[i]]));

			}
		}

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
		sprintf(line, "    %18s%18s\tStage%2d\n", sgpf.timeCols[t].c_str(), sgpf.timeRows[t].c_str(), t); tFile << line;
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
	int idx = 0;
	for (int c = 0; c < (int) sgpf.stocRows[0].size() ; c++ ) {
		if ( sgpf.stocCols[0][c] == "RHS" ) {
			sprintf(line, "    U%02d%14s%14.1f%14.1f\n", idx++, "LAGGED", 0.0, 0.02*vals);
			sFile << line;
		}
	}
	for (int c = 0; c < (int) sgpf.stocRows[0].size() ; c++ ) {
		if ( sgpf.stocRows[0][c] == "obj" ) {
			sprintf(line, "    U%02d%14s%14.1f%14.1f\n", idx++, "LAGGED", 0.0, 1.0);
			sFile << line;
		}
	}

	/* Establish the linear relationship between random variables of the stochastic process. */
	sprintf(line, "BLOCKS     LINTR\n"); sFile << line;

	/* Random variables in the model, constant terms based on past values. */
	for ( int t = 1; t < sgpf.numPeriods; t++ ) {
		sprintf(line, " BL PERIOD%d STAGE1\n", t); sFile << line;
		/* Root stage columns and rows do not have uncertain elements, hence index-0 for stocCols/stocRows corresponds
		 * to second stage. */
		for ( int i = 0; i < data.numMaturities; i++ ) { /* Computing the initial volume to be used for volume change brownian motion */
			vals += data.initVol[i];
		}
		for (int c = 0; c < (int) sgpf.stocRows[t-1].size(); c++ ) {
			if ( sgpf.stocCols[t-1][c] == "RHS" ) {
				sprintf(line, "    %14s%14s%14f\n", sgpf.stocCols[t-1][c].c_str(), sgpf.stocRows[t-1][c].c_str(), 0.001*vals);
				sFile << line;
			}
		}
		for (int c = 0; c < (int) sgpf.stocRows[t-1].size(); c++ ) {
			if ( sgpf.stocRows[t-1][c] == "obj" ) {
				if ( t == 1 )
					sprintf(line, "    %14s%14s%14f\n", sgpf.stocCols[t-1][c].c_str(), sgpf.stocRows[t-1][c].c_str(),
							data.initRet[c]);
				else
					sprintf(line, "    %14s%14s%14f\n", sgpf.stocCols[t-1][c].c_str(), sgpf.stocRows[t-1][c].c_str(),
							0.0);
				sFile << line;
			}
		}
	}

	idx = 0;
	/* Linear transformation matrix */
	for (int c = 0; c < (int) sgpf.stocRows[0].size(); c++ ) {
		if ( sgpf.stocCols[0][c] == "RHS" ) {
			sprintf(line, " RV U%02d%11s%14s\n", idx, "LAGGED", "LAG00"); sFile << line;
			sprintf(line, "    %14s%14s%14f\n", "RHS", sgpf.stocRows[0][c].c_str(), 1.0); sFile << line;
		}
	}
	for (int c = 0; c < (int) sgpf.stocCols[0].size(); c++ ) {
		if ( sgpf.stocRows[0][c] == sgpf.objName ) {
			sprintf(line, " RV U%02d%11s%14s\n", idx, "LAGGED", "LAG00"); sFile << line;
			sprintf(line, "    %14s%14s%14f\n", sgpf.stocCols[0][c].c_str(), sgpf.objName.c_str(), 1.0); sFile << line;
		}
	}

	for ( int t = 2 ; t < 3; t++ ) {
		for (int c = 0; c < (int) sgpf.stocRows[0].size(); c++ ) {
			if ( sgpf.stocCols[t-2][c] == "RHS" ) {
				sprintf(line, " HV %14s%14s%14s\n", "RHS", sgpf.stocRows[t-2][c].c_str(), "LAG01"); sFile << line;
				sprintf(line, "    %14s%14s%14f\n", "RHS", sgpf.stocRows[t-1][c].c_str(), 1.0); sFile << line;
			}
		}
		for (int c = 0; c < (int) sgpf.stocCols[0].size(); c++ ) {
			if ( sgpf.stocRows[t-2][c] == sgpf.objName ) {
				sprintf(line, " HV %14s%14s%14s\n", sgpf.stocCols[t-2][c].c_str(), sgpf.objName.c_str(), "LAG01"); sFile << line;
				sprintf(line, "    %14s%14s%14f\n", sgpf.stocCols[t-1][c].c_str(), sgpf.objName.c_str(), 1.0); sFile << line;
			}
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

	data.delta = 0.05;

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
			while ( t == 0 && ((data.retBorrow[t][n] = number) < CPP_TOL) );

			if ( t != 0 )
				while ( (data.retBorrow[t][n] = data.retBorrow[t-1][n] + number) < CPP_TOL )
					number = distribution(generator);
		}
	}

	/* Initial volume */
	data.initVol[1] = 200; data.initVol[2] = 120; data.initVol[3] = 80;

	return;
}//END defineData()
