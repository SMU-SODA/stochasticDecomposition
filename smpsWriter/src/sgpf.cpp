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

ILOSTLBEGIN

int numT, numMaturities, numAdq;
vector<double> ret, initVol;
vector<int> stdMaturities;
vector<vector<double>> retBorrow, retLend;

int createSGPFInstance() {

	/* Generate cor file */
	createSGPFcor();

	return 0;
}//END createSGPFInstance()

int createSGPFcor() {

	IloEnv env;
	char elemName[NAMESIZE];

	try {
		IloModel model(env);

		/* Setup model parameters */
		defineSGPFData(env);

		/* Decision variables: declaration */
		IloArray<IloNumVarArray> borrow(env, numT);
		IloArray<IloNumVarArray> lend(env, numT);
		IloArray<IloNumVarArray> vol(env, numT);
		IloNumVarArray totalVol(env, numT, 0, IloInfinity);

		/* Load the decision variables onto the solver */
		for ( int t = 0; t < numT; t++ ) {
			borrow[t] = IloNumVarArray (env, numMaturities, 0, IloInfinity);
			lend[t] = IloNumVarArray (env, numMaturities, 0, IloInfinity);
			vol[t] = IloNumVarArray (env, numMaturities, 0, IloInfinity);

			for ( int i = 0; i < numMaturities; i++ ) {
				/* Volume of maturity borrowed */
				sprintf(elemName, "borrow[%d][%d]", t, i);
				borrow[t][i].setName(elemName); model.add(borrow[t][i]);

				/* Volume of maturity lent */
				sprintf(elemName, "lend[%d][%d]", t, i);
				lend[t][i].setName(elemName); model.add(lend[t][i]);

				/* State/volume of the maturity. */
				sprintf(elemName, "vol[%d][%d]", t, i);
				vol[t][i].setName(elemName); model.add(vol[t][i]);
			}

			/* State of the portfolio */
			sprintf(elemName, "totalVol[%d]", t);
			totalVol[t].setName(elemName); model.add(totalVol[t]);
		}

		/* Objective function */
		IloObjective obj = IloMinimize(env, 0.0);
		IloExpr totalCost(env);
		for ( int t = 0; t < numT; t++) {
			for ( int i = 0; i < (int) stdMaturities.size(); i++ ) {
				totalCost += (retBorrow[t][stdMaturities[i]]*borrow[t][stdMaturities[i]] - retLend[t][stdMaturities[i]]*lend[t][i]);
			}
			totalCost += ret[t]*totalVol[t];
		}
		obj.setExpr(totalCost);
		model.add(obj);
		totalCost.end();

		/* Constraints */
		for ( int t = 0; t < numT; t++ ) {
			/* a. State dynamics for a standard maturity */
			for ( int i = 0; i < numMaturities; i++ ) {
				IloExpr expr (env);
				sprintf(elemName, "dynamicsStd[%d][%d]", t, i);

				auto it = find(stdMaturities.begin(), stdMaturities.end(), i);
				if ( it !=  stdMaturities.end() ) {
					if ( t != 0 ) {
						/* The bonds are one year closer to maturity, therefore a _(i+1)_ bond in previous time is now a _i_ time periods
						 * away from maturity. */
						expr = vol[t][i] - vol[t-1][i+1] - borrow[t][i] + lend[t][i];
						IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
					}
					else {
						expr = vol[t][i] - borrow[t][i] + lend[t][i];
						IloConstraint c(expr == initVol[i+1]); c.setName(elemName); model.add(c);
					}
				}
			}

			/* b. State dynamics for a non-standard maturity */
			for ( int i = 0; i < numMaturities; i++ ) {
				IloExpr expr (env);
				sprintf(elemName, "dynamicsNonStd[%d][%d]", t, i);

				vector<int>::iterator it = find(stdMaturities.begin(), stdMaturities.end(), i);
				if ( it ==  stdMaturities.end() ) {
					if ( t != 0) {
						expr = vol[t][i] - vol[t-1][i];
						IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
					}
					else {
						expr = vol[t][i];
						IloConstraint c(expr == initVol[i]); c.setName(elemName); model.add(c);
					}
				}
			}

			/* c. State/volume of the portfolio */
			{
				IloExpr expr (env);
				sprintf(elemName, "pfState[%d]", t);

				expr = totalVol[t];
				for ( int i = 0; i < numMaturities; i++ ) {
					expr -= vol[t][i];
				}
				IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
			}

			/* d. Change in portfolio volume change */
			{
				IloExpr expr (env);
				sprintf(elemName, "volChange[%d]", t);

				double initTotalVol = 0.0;
				expr = totalVol[t];
				if ( t == 0 ) {
					for ( int i = 0; i < numMaturities; i++ ) {
						initTotalVol += initVol[i];
					}
				}
				else {
					expr -= totalVol[t-1];
				}
				IloConstraint c(expr == initTotalVol); c.setName(elemName); model.add(c);
			}

			/* e. Adequacy inequality */
			{
				IloExpr expr (env);
				sprintf(elemName, "adequacy[%d]", t);

				for ( int i = 0; i < numMaturities; i++ ) {
					vector<int>::iterator it = find(stdMaturities.begin(), stdMaturities.end(), i);
					if ( it ==  stdMaturities.end() ) {
						expr += borrow[t][i];
					}
				}
				double adqTotalVol = 0.0;

				if ( t != 0 ) {
					for ( int i = 0; i < numAdq; i++ ) {
						expr -= vol[t-1][i];
					}
				}
				else {
					for ( int i = 0; i < numAdq; i++ ) {
						adqTotalVol += initVol[i];
					}
				}
				IloConstraint c(expr <= adqTotalVol); c.setName(elemName); model.add(c);
			}
		}

		IloCplex cplex(model);
		cplex.exportModel("sgpf.lp");
		cplex.solve();
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	env.end();

	return 0;
}//END createSGPF()

void defineSGPFData(IloEnv env) {

	numT = 3;			/* Number of time periods */
	numMaturities = 60;	/* Total number of maturities */
	numAdq = 3;			/* Number of maturities used to ensure adequacy (M in the paper) */

	stdMaturities = {0,1,2,3,4,5,6,7,8}; /* list of standard maturities */

	ret = {0.005, 0.004978716, 0.004957761, 0.00493713, 0.004916816, 0.004896814};	/* Return on the total volume */

	/* Return rate for borrowing */
	retBorrow = {{0.004281696,0.008497303,0.01264803, 0.024722558,0.024190509,0.023507712,0.022532289,0.021556865,0.020581441},
			{0.005337856,0.010332067,0.01498735, 0.022560369,0.018126625,0.018065661,0.017760841,0.017456021,0.017151201},
			{0.006385355,0.012151986,0.017308013,0.019614885,0.012875594,0.013233249,0.01339582,0.01355839,0.013720961},
			{0.007424393,0.013957405,0.019610451,0.015886106,0.008437415,0.009010477,0.009437225,0.009863973,0.010290721},
			{0.008455166,0.01574866 ,0.014655003,0.011374032,0.004812091,0.005397345,0.005885057,0.006372769,0.00686048},
			{0.009477867,0.008798027,0.008118186,0.006078663,0.001999619,0.002393852,0.002739315,0.003084778,0.00343024}};

	/* Return rate on lending */
	retLend = {{0.004231909,0.008398147,0.012499916,0.024429931,0.023897882,0.022922458,0.021947034,0.020971611,0.019996187},
			{0.005288279,0.010233324,0.014839846,0.022316513,0.017882769,0.017577949,0.017273129,0.016968309,0.016663489},
			{0.006335984,0.01205365,0.017161109,0.0194198,0.012680509,0.012843079,0.01300565,0.013168221,0.013330791},
			{0.007375225,0.013859469,0.019464138,0.015739792,0.008291102,0.00871785,0.009144598,0.009571346,0.009998093},
			{0.008406198,0.015651118,0.014557461,0.01127649,0.004714548,0.00520226,0.005689972,0.006177684,0.006665396},
			{0.009429096,0.008749256,0.008069415,0.006029892,0.001950848,0.00229631,0.002641773,0.002987235,0.003332698}};

	/* Initial volume. There are 61 entries here. The first entry matures before the model horizon begins, but is retained in the
	 * data set to maintain consistency. */
	initVol = {200000, 120000, 80000, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

}//END defineData()
