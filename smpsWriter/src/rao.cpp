
#include "writer.hpp"
#include <stdlib.h>

class RaoData {

public:
	int numPeriods;
	int numProducts;
	bool fixed;
	vector<vector<double>> meanDemand;
	vector<vector<double>>  stdDemand;
	vector<vector<double>> overageCost, shortageCost, productionCost;
	vector<vector<vector<double>>> subCost;
};

void defineRaoData(SMPSmodel& Rao1, RaoData& data);
int createRaoCor(SMPSmodel& Rao1, RaoData& data, char probName[128]);
int createRaoTim(SMPSmodel Rao, char probName[128]);
int createRaoStoc(SMPSmodel Rao, RaoData& data, char probName[128]);

int createRaoInstance(string inputDir, string outputDir, int argc, char* argv[]) {
	SMPSmodel Rao;
	RaoData data;
	char srcFile[128], destnFile[128], probName[64];

	/* Setup model parameters */
	for (int i = 1; (i < argc); i++) {
		if (argv[i][0] == '-') {
			switch ((argv[i])[1]) {
			case 't': data.numPeriods = atoi(argv[++i]); break;
			case 'n': data.numProducts= atoi(argv[++i]); break;
			case 'f': data.fixed= atoi(argv[++i]); break;
			default: ++i; continue; 
			}
		}
		else {
			printf("Input options (%s) must begin with a '-'. Use '-?' for help.\n", argv[i]); exit(0);
		}
	}

	sprintf(probName, "raoP%02dT%02d", data.numProducts, data.numPeriods);

	defineRaoData(Rao, data);

	/* Generate cor file */
	createRaoCor(Rao, data, probName);
	createRaoTim(Rao, probName);
	createRaoStoc(Rao, data, probName);

	/* Move the files to the spInput folder */
	sprintf(destnFile, "%s.cor", probName);
	sprintf(srcFile, "%s.mps", probName);
	rename(srcFile, destnFile);

	sprintf(destnFile, "mkdir %s%s\\", outputDir.c_str(), probName);
	system(destnFile);

	sprintf(destnFile, "%s%s\\", outputDir.c_str(), probName);
#ifdef _WIN64
	sprintf(srcFile, "move %s.* %s", probName, destnFile);
#else
	sprintf(srcFile, "mv %s.* %s", probName, destnFile);
#endif
	system(srcFile);

	printf("Successfully wrote the SMPS files for %s to the output directory '%s'.\n", probName, outputDir.c_str());

	return 0;
}//END createRaoInstance()

int createRaoCor(SMPSmodel& Rao, RaoData& data, char probName[128]) {
	char elemName[NAMESIZE];

	try {
		IloEnv   env;
		IloModel model(env, probName);

		/**************** Decision variables *****************/
		IloArray<IloNumVarArray> prodQuantity(env, data.numPeriods);
		IloArray<IloNumVarArray> excess(env, data.numPeriods);
		IloArray<IloNumVarArray> shortage(env, data.numPeriods);
		IloArray<IloArray<IloNumVarArray>> substitution(env, data.numPeriods);

		/* Loop over time periods to make sure that columns are time stage ordered */
		/* Fisrt stage decisions*/
		for (int t = 0; t < data.numPeriods; t++) {
			/* Generating the production quantity decision variables for a given time period*/
			prodQuantity[t] = IloNumVarArray(env, data.numProducts, 0, 3000, ILOFLOAT);
			for (int i = 0; i < data.numProducts; i++) {
				sprintf(elemName, "prodQuantity[%d][%d]", t, i);
				prodQuantity[t][i].setName(elemName); model.add(prodQuantity[t][i]);
				if (i == 0 && t == 0) {
					Rao.timeCols.push_back(elemName);
				}
			}
		}
		/* Second stage decisions*/
		for (int t = 0; t < data.numPeriods; t++) {
			/* Generating the substitution decision variables for a given time period*/
			substitution[t] = IloArray<IloNumVarArray>(env, data.numProducts);
			/*Looping through the substituting product */
			for (int i = 0; i < data.numProducts; i++) {
				substitution[t][i] = IloNumVarArray(env, data.numProducts, 0, IloInfinity, ILOFLOAT);
				/* Looping through the substituted product*/
				for (int j = i; j < data.numProducts; j++) {
					sprintf(elemName, "substitution[%d][%d][%d]", t, i, j);
					substitution[t][i][j].setName(elemName); model.add(substitution[t][i][j]);
					if (i == 0 && j == 0 && t == 0) {
						Rao.timeCols.push_back(elemName);
					}

				}
			}

			/* Generating the excess decision variables for a given time period*/
			excess[t] = IloNumVarArray(env, data.numProducts, 0, IloInfinity, ILOFLOAT);
			for (int i = 0; i < data.numProducts; i++) {
				sprintf(elemName, "excess[%d][%d]", t, i);
				excess[t][i].setName(elemName); model.add(excess[t][i]);
			}
			/* Generating the shortage decision variables for a given time period*/
			shortage[t] = IloNumVarArray(env, data.numProducts, 0, IloInfinity, ILOFLOAT);
			for (int i = 0; i < data.numProducts; i++) {
				sprintf(elemName, "shortage[%d][%d]", t, i);
				shortage[t][i].setName(elemName); model.add(shortage[t][i]);
			}
		}

		/**************** Objective function *****************/

		IloObjective obj = IloMinimize(env, 0.0);
		IloExpr totalCost(env);

		for (int t = 0; t < data.numPeriods; t++) {
			for (int i = 0; i < data.numProducts; i++) {
				/* Adding production costs, shortage costs and overage costs to the objective function*/
				totalCost += data.productionCost[t][i] * prodQuantity[t][i] + data.shortageCost[t][i] * shortage[t][i] + data.overageCost[t][i] * excess[t][i];
				/* Adding substitution costs to the objective function*/
				for (int j = i; j < data.numProducts; j++) {
					totalCost += data.subCost[t][i][j] * substitution[t][i][j];
				}
			}
		}
		Rao.objName = "obj";
		obj.setExpr(totalCost);
		obj.setName("obj");
		model.add(obj); totalCost.end();
		Rao.timeRows.push_back("obj");

		/**************** Constraints *****************/


		/* Loop over time periods to make sure constraints are time stage ordered */
		for (int t = 0; t < data.numPeriods; t++) {

			/* Supply balance constraint for each product at a given time period*/
			for (int i = 0; i < data.numProducts; i++) {
				IloExpr expr(env);
				sprintf(elemName, "SupplyBalance[%d][%d]", t, i);
				for (int j = i; j < data.numProducts; j++) {
					expr += substitution[t][i][j];
				}
				expr += excess[t][i];
				expr += -1 * prodQuantity[t][i];
				if (t > 0) {
					expr += -1 * excess[t - 1][i];
				}
				Rao.timeRows.push_back(elemName);
				IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
			}

			/* Demand balance constraint for each product at a given time period*/
			for (int j = 0; j < data.numProducts; j++) {
				IloExpr expr(env);
				sprintf(elemName, "DemandBalance[%d][%d]", t, j);
				for (int i = 0; i <= j; i++) {
					expr += substitution[t][i][j];
				}
				expr += shortage[t][j];
				if (t > 0) {
					expr += -1 * shortage[t - 1][j];
				}
				Rao.stocRows[0].push_back(elemName);
				Rao.stocCols[0].push_back(" RHS");
				IloConstraint c(expr == data.meanDemand[t][j]); c.setName(elemName); model.add(c);
			}
		}
		IloCplex cplex(model);
		sprintf(elemName, "%s.mps", probName);
		cplex.exportModel(elemName);
		sprintf(elemName, "%s.lp", probName);
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

}//END createRaoCor()

int createRaoTim(SMPSmodel Rao, char probName[128]) {
	ofstream tFile;
	char fName[NAMESIZE], line[BLOCKSIZE];

	sprintf(fName, "%s.tim", probName);
	tFile.open(fName);
	tFile << "TIME   " << probName << endl;

	tFile << "PERIODS" << endl;
	for (int t = 0; t < Rao.numStages; t++) {  //??????????????????????????????????????????[t]
		sprintf(line, "    %25s%25s\tTIME%02d\n", Rao.timeCols[t].c_str(), Rao.timeRows[t].c_str(), t); tFile << line;
	}
	tFile << "ENDATA" << endl;
	tFile.close();

	return 0;
}//END createRaoTim()

int createRaoStoc(SMPSmodel Rao, RaoData& data, char probName[128]) {
	ofstream sFile;
	char fName[NAMESIZE], line[BLOCKSIZE];

	sprintf(fName, "%s.sto", probName); sFile.open(fName);

	sprintf(line, "%-14s%-14s\n", "STOCH", probName);
	sFile << line;
	sFile << "INDEP          NORMAL" << endl;
	int c = 0;
	for (int t = 0; t < data.numPeriods; t++) {
		for (int i = 0; i < data.numProducts; i++) {
			if (Rao.stocCols[0][c] == " RHS") {
				sprintf(line, "%s%25s%25f%25f\n", Rao.stocCols[0][c].c_str(), Rao.stocRows[0][c].c_str(), data.meanDemand[t][i], data.stdDemand[t][i]);

				sFile << line;
				c++;
			}
		}
	}
	sFile << "ENDATA" << endl;
	sFile.close();

	return 0;
}//END createRaoStoc()

void defineRaoData(SMPSmodel& Rao, RaoData& data) {
	data.subCost = vector<vector<vector<double>>>(data.numPeriods, vector<vector<double>>(data.numProducts, vector<double>(data.numProducts)));
	data.overageCost = vector<vector<double>>(data.numPeriods, vector<double>(data.numProducts));
	data.productionCost = vector<vector<double>>(data.numPeriods, vector<double>(data.numProducts));
	data.shortageCost = vector<vector<double>>(data.numPeriods, vector<double>(data.numProducts));
	data.meanDemand = vector<vector<double>>(data.numPeriods, vector<double>(data.numProducts));
	data.stdDemand = vector<vector<double>>(data.numPeriods, vector<double>(data.numProducts));

	
	/* Loop over time periods */
	for (int t = 0; t < data.numPeriods; t++) {

		/* The production cost is calculated using the Rao's paper fomulation :
	        productionCost[i]= (n-1-i)* \eta +1; where n is the number of products and \eta can be 0.1; 0.2; 0.5  */
		
		for (int i = 0; i < data.numProducts; i++) {
			data.productionCost[t][i] = double((data.numProducts - 1 - i) * 0.2 + 1.0);
		}


		/* The substitution cost is definded only for the cases where i<j
		   and is calculated using the Rao's paper fomulation:
		   (the different between their corresponding production cost)* T ; where T = 0.1; 0.5; 1  */

		for (int i = 0; i < data.numProducts; i++) {
			for (int j = i+1 ; j < data.numProducts; j++) {
				data.subCost[t][i][j] = double((j - i) * 0.1);
			}
		}


		/* The overage cost is the holding cost for all time periods except the last time period
			which is the difference between holding cost and salvage value, the holding cost assumed to be
			20 percent of the production cost and for the overage cost the Rao's paper formulation is used
			where gi-si = -65,-70, 0, 15 percent of ci */

		for (int i = 0; i < data.numProducts; i++) {
			if (t == data.numPeriods - 1) {
				data.overageCost[t][i] = -0.65 * data.productionCost[t][i];
			}
			else {
				data.overageCost[t][i] = 0.2 * data.productionCost[t][i];
			}

		}

		/* The shortage value is calculated considering a critical rati; where (pi-ci)/(pi+(gi-si))= 0.5, 0.7, 0.8, 0.85, 0.9, 0.95 , pi: shortage cost */

		for (int i = 0; i < data.numProducts; i++) {
			data.shortageCost[t][i] = (data.productionCost[t][i] + data.overageCost[data.numPeriods - 1][i] * 0.85) / (1 - 0.85);
			data.stdDemand[t][i] = int((double)rand() / (RAND_MAX + 1) * (40 - 20)+ 20);
			data.meanDemand[t][i] = int((double)rand() / (RAND_MAX + 1) * (140 - 70) + 70);
		}



		
	}
	/* prepare to store the columns and rows with random elements in them. */
	Rao.stocCols = vector<vector<string>>(1, vector<string>());
	Rao.stocRows = vector<vector<string>>(1, vector<string>());
}




