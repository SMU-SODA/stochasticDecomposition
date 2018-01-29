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

void defineData(IloEnv env);

ILOSTLBEGIN

int numT, numMaturities, numAdq;
vector<double> ret;
vector<int> stdMaturities;
vector<vector<double>> retBorrow, retLend, initVol;

int createSGPF() {
	IloEnv env;
	char elemName[NAMESIZE];

	try {
		IloModel model(env);

		/* Setup model parameters */
		defineData(env);

		/* Decision variables: declaration */
		IloArray<IloNumVarArray> borrow(env, numT);
		IloArray<IloNumVarArray> lend(env, numT);
		IloArray<IloNumVarArray> vol(env, numT);
		IloNumVarArray totalVol(env, numT);

		/* Load the decision variables onto the solver */
		for ( int t = 0; t < numT; t++ ) {
			borrow[t] = IloNumVarArray > (env, numMaturities, 0, IloInfinity);
			lend[t] = IloNumVarArray > (env, numMaturities, 0, IloInfinity);
			vol[t] = IloNumVarArray > (env, numMaturities);

			for ( int i = 0; i < numMaturities; i++ ) {
				/* Volume of maturity borrowed */
				sprintf(elemName, "borrow[%d][%d]", t, i);
				borrow[t][i].setName(elemName); model.add(borrow[t][i]);

				/* Volume of maturity lent */
				sprintf(elemName, "lend[%d][%d]", t, i);
				lend[t][i].setName(elemName); model.add(lend[t][i]);

				/* State/volume of the maturity. */
				sprintf(elemName, "vol[%d]", t);
				vol[t][i].setName(elemName); model.add(vol[i][t]);
			}

			/* State of the portfolio */
			sprintf(elemName, "totalVol[%d]", t);
			totalVol[t].setName(elemName); model.add(totalVol[t]);
		}

		/* Objective function */
		IloObjective obj;
		IloExpr totalCost(env);
		for ( int t = 0; t < numT; t++) {
			for ( int i = 0; i < stdMaturities.size(); i++ ) {
				totalCost += (retLend[t][stdMaturities[i]]*lend[t][i] - retBorrow[t][stdMaturities[i]]*borrow[t][stdMaturities[i]]);
			}
			totalCost += ret[t]*vol[t];
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
						expr = vol[t][i] - vol[t-1][i] - borrow[t][i] + lend[t][i];
						IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
					}
					else {
						expr = vol[t][i] - borrow[t][i] + lend[t][i];
						IloConstraint c(expr == initVol[i]); c.setName(elemName); model.add(c);
					}
				}
			}

			/* b. State dynamics for a non-standard maturity */
			for ( int i = 0; i < numMaturities; i++ ) {
				IloExpr expr (env);
				sprintf(elemName, "dynamicsNonStd[%d][%d]", t, i);

				vector<int>::iterator it = find(stdMaturities.begin(), stdMaturities.end(), i);
				if ( it ==  stdMaturities.end() ) {
					if ( t == 0) {
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

				if ( t == 0 ) {
					expr = totalVol[t];
				}
				else {
					expr = totalVol[t] - totalVol[t-1];
				}
				IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
			}

			/* e. Adequacy inequality */
			{
				IloExpr expr (env);
				sprintf(elemName, "adequacy[%d]", t);

				if ( t == 0 ) {
					for ( int i = 0; i < numMaturities; i++ ) {
						vector<int>::iterator it = find(stdMaturities.begin(), stdMaturities.end(), i);
						if ( it ==  stdMaturities.end() ) {
							expr += borrow[t][i];
						}
					}

					for ( int i = 0; i < numAdq; i++ ) {
						expr -= vol[t][i];
					}
				}
				else {
					expr = totalVol[t] - totalVol[t-1];
				}
				IloConstraint c(expr == 0); c.setName(elemName); model.add(c);

			}
		}
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

void defineData(IloEnv env) {

	numT = 3;			/* Number of time periods */
	numMaturities = 60;	/* Total number of maturities */
	numAdq = 3;			/* Number of maturities used to ensure adaquacy (M in the paper) */

	stdMaturities = {};

	vector<double> ret;
	vector<int> ;
	vector<vector<double>> retBorrow, retLend, initVol;


}//END defineData()
