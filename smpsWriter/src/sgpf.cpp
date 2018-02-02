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
};

void defineSGPFdata(SMPSmodel &sgpf, SGPFdata &data);
int createSGPFcor(SMPSmodel &sgpf, SGPFdata data);
int createSGPFtim(SMPSmodel sgpf);
int createSGPFstoc (SMPSmodel sgpf, SGPFdata data);

int createSGPFInstance() {
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

	sprintf(destnFile, "mkdir %s/sgpf%dy%d/", sgpf.dir.c_str(), sgpf.numPeriods, sgpf.numStages);
	system(destnFile);

	sprintf(srcFile, "mv sgpf%dy%d.* %s", sgpf.numPeriods, sgpf.numStages, destnFile);
	system(srcFile);

	return 0;
}//END createSGPFInstance()

int createSGPFcor(SMPSmodel &sgpf, SGPFdata data) {
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

			for ( int i = 0; i < data.numMaturities; i++ ) {
				vector<int>::iterator it = find(data.stdMaturities.begin(), data.stdMaturities.end(), i);
				if ( it !=  data.stdMaturities.end() ) {
					/* Volume of maturity borrowed */
					sprintf(elemName, "borrow[%d][%d]", t, i);
					borrow[t][i].setName(elemName); model.add(borrow[t][i]);
					if ( i == 0 )
						sgpf.timCols.push_back(elemName);
					if ( t != 0 )
						sgpf.stocCols.push_back(elemName);

					/* Volume of maturity lent */
					sprintf(elemName, "lend[%d][%d]", t, i);
					lend[t][i].setName(elemName); model.add(lend[t][i]);
					if ( t != 0 )
						sgpf.stocCols.push_back(elemName);
				}

				/* State/volume of the maturity. */
				sprintf(elemName, "volume[%d][%d]", t, i);
				volume[t][i].setName(elemName); model.add(volume[t][i]);
			}

			/* State of the portfolio */
			sprintf(elemName, "totalVolume[%d]", t);
			totalVolume[t].setName(elemName); model.add(totalVolume[t]);
			if ( t != 0 )
				sgpf.stocCols.push_back(elemName);
		}

		/* Objective function */
		sgpf.objName = "netReturn";
		IloObjective obj = IloMinimize(env, 0.0, "netReturn");
		IloExpr netReturn(env);
		for ( int t = 0; t < sgpf.numPeriods; t++) {
			for ( int i = 0; i < data.numStd; i++ ) {
				netReturn += (data.retBorrow[t][data.stdMaturities[i]]*borrow[t][data.stdMaturities[i]] - data.retLend[t][data.stdMaturities[i]]*lend[t][i]);
			}
			netReturn += data.ret[t]*totalVolume[t];
		}
		obj.setExpr(netReturn); model.add(obj); netReturn.end();

		/* Constraints */
		for ( int t = 0; t < sgpf.numPeriods; t++ ) {
			/* a. State dynamics for a standard maturity */
			for ( int i = 0; i < data.numMaturities; i++ ) {
				IloExpr expr (env);
				sprintf(elemName, "dynamicsStd[%d][%d]", t, i);

				vector<int>::iterator it = find(data.stdMaturities.begin(), data.stdMaturities.end(), i);
				if ( it !=  data.stdMaturities.end() ) {
					if ( t != 0 ) {
						/* The bonds are one year closer to maturity, therefore a _(i+1)_ bond in previous time is now a _i_ time periods
						 * away from maturity. */
						expr = volume[t][i] - volume[t-1][i+1] - borrow[t][i] + lend[t][i];
						IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
					}
					else {
						expr = volume[t][i] - borrow[t][i] + lend[t][i];
						IloConstraint c(expr == data.initVol[i+1]); c.setName(elemName); model.add(c);
					}
				}
				if ( i == 0 )
					sgpf.timRows.push_back(elemName);
			}

			/* b. State dynamics for a non-standard maturity */
			for ( int i = 0; i < data.numMaturities; i++ ) {
				IloExpr expr (env);
				sprintf(elemName, "dynamicsNonStd[%d][%d]", t, i);

				vector<int>::iterator it = find(data.stdMaturities.begin(), data.stdMaturities.end(), i);
				if ( it ==  data.stdMaturities.end() ) {
					if ( t != 0) {
						expr = volume[t][i] - volume[t-1][i];
						IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
					}
					else {
						expr = volume[t][i];
						IloConstraint c(expr == data.initVol[i]); c.setName(elemName); model.add(c);
					}
				}
			}

			/* c. State/volume of the portfolio */
			{
				IloExpr expr (env);
				sprintf(elemName, "pfState[%d]", t);

				expr = totalVolume[t];
				for ( int i = 0; i < data.numMaturities; i++ ) {
					expr -= volume[t][i];
				}
				IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
			}

			/* d. Change in portfolio volume change */
			{
				IloExpr expr (env);
				sprintf(elemName, "volChange[%d]", t);

				sgpf.stocRows.push_back(elemName);

				double initTotalVol = 0.0;
				expr = totalVolume[t];
				if ( t == 0 ) {
					for ( int i = 0; i < data.numMaturities; i++ ) {
						initTotalVol += data.initVol[i];
					}
				}
				else {
					expr -= totalVolume[t-1];
				}
				IloConstraint c(expr == initTotalVol); c.setName(elemName); model.add(c);
			}

			/* e. Adequacy inequality */
			{
				IloExpr expr (env);
				sprintf(elemName, "adequacy[%d]", t);

				sgpf.stocRows.push_back(elemName);

				for ( int i = 0; i < data.numMaturities; i++ ) {
					vector<int>::iterator it = find(data.stdMaturities.begin(), data.stdMaturities.end(), i);
					if ( it !=  data.stdMaturities.end() ) {
						expr += borrow[t][data.stdMaturities[i]];
					}
				}
				double adqTotalVol = 0.0;

				if ( t != 0 ) {
					for ( int i = 0; i < data.numAdq; i++ ) {
						expr -= volume[t-1][i];
					}
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
	char fName[NAMESIZE];

	sprintf(fName, "sgpf%dy%d.tim", sgpf.numPeriods, sgpf.numStages);
	tFile.open(fName);
	tFile << "TIME   sgpf" << sgpf.numPeriods << "y" << sgpf.numStages << endl;

	tFile << "PERIODS" << endl;
	for (int t = 0; t < sgpf.numStages; t++ )
		tFile << setw(4) << setw(14) << sgpf.timCols[t] << "\t" << sgpf.timRows[t] << endl;
	tFile << "ENDATA" << endl;
	tFile.close();

	return 0;
}//END createSGPFtim()

int createSGPFstoc (SMPSmodel sgpf, SGPFdata data) {
	ofstream sFile;
	char fName[NAMESIZE];
	double vals;

	vals = 0.0;
	for ( int i = 0; i < data.numMaturities; i++ ) {
		vals += data.initVol[i];
	}

	sprintf(fName, "sgpf%dy%d.sto", sgpf.numPeriods, sgpf.numStages);
	sFile.open(fName);
	sprintf(fName, "sgpf%dy%d", sgpf.numPeriods, sgpf.numStages);
	sFile << "STOCH"  << setw(14) << fName << endl;
	sFile << "INDEP          WEINER" << endl;
	for (int c = 0; c < (int) sgpf.stocCols.size(); c++ ) {
		sFile << "    " << setw(14) << sgpf.stocCols[c] << setw(14) << sgpf.objName << "\t0\t1" << endl;
	}
	for (int r = 0; r < (int) sgpf.stocRows.size(); r++ ) {
		sFile << "    " << setw(14) << "RHS" << setw(14) << sgpf.stocRows[r] << "\t" << 0.001*vals << "\t" << 0.00005*vals << endl;
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
		sgpf.numPeriods = 2;
		sgpf.numStages = 2;
	}

	data.numMaturities = 60;	/* Total number of maturities */
	data.numAdq = 3;			/* Number of maturities used to ensure adequacy (M in the paper) */
	data.stdMaturities = {0,1,2,3,4,5,6,7,8}; /* list of standard maturities */
	data.numStd = (int) data.stdMaturities.size();

	/* Return rate for return, borrowing, lending */
	data.ret = vector<double> (sgpf.numPeriods);
	data.retBorrow = vector< vector<double>> (sgpf.numPeriods, vector<double> (data.numStd) );
	data.retLend = vector< vector<double>> (sgpf.numPeriods, vector<double> (data.numStd) );
	data.initVol = vector<double> (data.numMaturities);

	default_random_engine generator{static_cast<long unsigned int>(3554548844580680)};
	normal_distribution <double> distribution(0.0,1.0);

	/* Returns */
	for ( int t = 0; t < sgpf.numPeriods; t++ ) {
		double number;
		do {
			number = distribution(generator);
		}
		while ( t == 0 && ((data.ret[t] = number) < TOLERANCE) );

		if ( t > 0 )
			while ( (data.ret[t] = data.ret[t-1] + number) < TOLERANCE )
				number = distribution(generator);
	}

	/* Borrowing rate */
	for ( int n = 0; n < data.numStd; n++ ) {
		for ( int t = 0; t < sgpf.numPeriods; t++ ) {
			double number;
			do {
				number = distribution(generator);
			}
			while ( t == 0 && ((data.retBorrow[t][data.stdMaturities[n]] = number) < TOLERANCE) );

			if ( t != 0 )
				while ( (data.retBorrow[t][data.stdMaturities[n]] = data.retBorrow[t-1][data.stdMaturities[n]] + number) < TOLERANCE )
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
			while ( t == 0 && ((data.retLend[t][data.stdMaturities[n]] = number) < TOLERANCE) );

			if ( t != 0 )
				while ( (data.retLend[t][data.stdMaturities[n]] = data.retBorrow[t-1][data.stdMaturities[n]] + number) < TOLERANCE )
					number = distribution(generator);
		}
	}

	/* Initial volume */
	data.initVol[0] = 200000; data.initVol[1] = 120000; data.initVol[2] = 80000;

	return;
}//END defineData()
