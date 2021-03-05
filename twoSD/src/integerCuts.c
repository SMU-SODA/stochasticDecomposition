///*
// * integerCuts.c
// *
// *  Created on: Feb 24, 2021
// *      Author: Harsha Gangammanavar
// * Institution: Southern Methodist University
// *  
// * Please send your comments or bug report to harsha (at) smu (dot) edu
// *
// */
//
//#include "twoSD.h"
//
//int addMIPCut2Pool(cellType *cell, oneCut *cut, int lenX, double lb, bool GMI);
//
// /*
//
// by siavash tabrizian Feb 2020
//
// subroutin for adding the MIP (GMI and MIR) cuts to their corresponding cutType structure
//
// */
//int addMIPCut2Pool(cellType *cell, oneCut *cut, int lenX, double lb, bool GMI) {
//
//	if (GMI) {
//		if (cell->GMIcuts->cnt <= cell->maxMIPCuts)
//		{
//			cell->GMIcuts->vals[cell->GMIcuts->cnt] = cut;
//			return cell->GMIcuts->cnt;
//			printf("%i \n", cell->GMIcuts->cnt);
//		}
//		else {
//			return 0;
//		}
//	}
//	else {
//		if (cell->MIRcuts->cnt <= cell->maxMIPCuts)
//		{
//			cell->MIRcuts->vals[cell->MIRcuts->cnt] = cut;
//			return cell->MIRcuts->cnt;
//		}
//		else {
//			return 0;
//		}
//	}
//
//}//END addCut()
//
// //----------------------------------------------------
// /*
//
// by siavash tabrizian Feb 20
//
// Form MIP cuts (GMI and MIR)
//
// */
//
//int formGMICut(probType **prob, cellType *cell, dVector Xvect, double lb) {
//	oneCut 	**cut;
//	int    	cutIdx;
//
//
//	/* (b) create a GMI cut */
//	clock_t tic = clock();
//	int cutnum_before = cell->GMIcuts->cnt;
//	cut = purefracGMICut(prob, cell, Xvect, lb);
//	int cutnum_after = cell->GMIcuts->cnt;
//	int newcuts = cutnum_after - cutnum_before;
//	if (cut == NULL) {
//		errMsg("algorithm", "formGMICut", "failed to create the affine minorant", 0);
//		return -1;
//	}
//	cell->time.argmaxIter += ((double)(clock() - tic)) / CLOCKS_PER_SEC;
//
//	for (int c = 0; c < newcuts; c++)
//	{
//		/* (c) add cut to the structure and master problem  */
//		if ((cutIdx = addMIPCut2Pool(cell, cut[c], prob[0]->num->cols, lb, true)) < 0) {
//			errMsg("algorithm", "formGMICut", "failed to add the new cut to cutsType structure", 0);
//			return -1;
//		}
//
//		if (addMIPCut2Master(cell->master, cut[c], cell->incumbX, prob[0]->num->cols, true)) {
//			errMsg("algorithm", "formGMICut", "failed to add the new cut to master problem", 0);
//			return -1;
//		}
//	}
//
//
//
//	cutIdx = cell->GMIcuts->cnt;
//
//	return cutIdx;
//}//END formCut()
//
//int formMIRCut(probType **prob, cellType *cell, dVector Xvect, double lb) {
//	oneCut 	**cut;
//	int    	cutIdx;
//
//
//	/* (b) create a MIR cut */
//	clock_t tic = clock();
//	int cutnum_before = cell->MIRcuts->cnt;
//	cut = pureMIRCut(prob, cell, Xvect, lb);
//	int cutnum_after = cell->MIRcuts->cnt;
//	int newcuts = cutnum_after - cutnum_before;
//	if (cut == NULL) {
//		errMsg("algorithm", "formGMICut", "failed to create the affine minorant", 0);
//		return -1;
//	}
//	cell->time.argmaxIter += ((double)(clock() - tic)) / CLOCKS_PER_SEC;
//
//	for (int c = 0; c < newcuts; c++)
//	{
//		/* (c) add cut to the structure and master problem  */
//		if ((cutIdx = addMIPCut2Pool(cell, cut[c], prob[0]->num->cols, lb, false)) < 0) {
//			errMsg("algorithm", "formGMICut", "failed to add the new cut to cutsType structure", 0);
//			return -1;
//		}
//
//		if (addMIPCut2Master(cell->master, cut[c], cell->incumbX, prob[0]->num->cols, false)) {
//			errMsg("algorithm", "formGMICut", "failed to add the new cut to master problem", 0);
//			return -1;
//		}
//	}
//
//
//
//	cutIdx = cell->MIRcuts->cnt;
//
//	return cutIdx;
//}//END formCut()
//
////
//// /*
////
//// form the pure integer MIR cuts
////
//// siavash tabrizian 30 Nov
////
//// */
////oneCut **pureMIRCut(probType **prob, cellType *cell, dVector Xvect, double lb) {
////	oneCut 	*cut;
////	oneCut 	**cutarr;
////	int 	k;
////
////	if (!(cutarr = (oneCut **) arr_alloc(cell->cuts->cnt + cell->cuts->cnt*cell->cuts->cnt / 2, oneCut *)))
////		errMsg("allocation", "formMIRcut", "cutarr", 0);
////
////
////	int cutnum = 0;
////	//// Creating MIR for benders cuts
////	for (k = 0; k < cell->cuts->cnt; k++) {
////
////		/* allocate memory to a new MIR cut */
////		cut = newCut(prob[0]->num->cols, 0, 0);
////
////		double bfloor = floor(cell->cuts->vals[k]->alpha);
////		double bceil = ceil(cell->cuts->vals[k]->alpha);
////		double f0 = cell->cuts->vals[k]->alpha - bfloor;
////		cut->beta[0] = 1;
////		cut->alpha = bceil * f0;
////
////		for (int v = 1; v < cell->master->mac; v++)
////		{
////
////			double afloor = floor(cell->cuts->vals[k]->beta[v]);
////			double aceil = floor(cell->cuts->vals[k]->beta[v]);
////			double fj = cell->cuts->vals[k]->beta[v] - afloor;
////			double beta1 = afloor * f0 + fj;
////			double beta2 = aceil * f0;
////
////			cut->beta[v] = minimum(beta1, beta2);
////
////		}
////
////		bool repeatFlag = false;
////		if (cell->MIRcuts->vals != NULL)
////		{
////			for (int pc = 0; pc < cell->MIRcuts->cnt; pc++)
////			{
////				if (cell->MIRcuts->vals[pc] != NULL && cut->alpha == cell->MIRcuts->vals[pc]->alpha &&
////					twoNorm(cut->beta, cell->MIRcuts->vals[pc]->beta, cell->master->mac) <= 0.001)
////				{
////					repeatFlag = true;
////					break;
////				}
////			}
////		}
////
////		if (repeatFlag == false)
////		{
////			cut->numSamples = cell->cuts->vals[k]->numSamples;
////			cell->MIRcuts->cnt++;
////			cutarr[cutnum] = cut;
////			cutnum++;
////		}
////
////	}
////
////
////#if defined(MIRSubbaddActive)
////	//// Creating subbadditive cuts from original MIRs
////	for (k = 0; k < cell->cuts->cnt - 1; k++)
////	{
////		for (int k2 = k + 1; k2 < cell->cuts->cnt; k2++)
////		{
////			/* allocate memory to a new MIR cut */
////			cut = newCut(prob[0]->num->cols, 0, 0);
////			//printf("here\n");
////			double alphfloor1 = floor(cell->cuts->vals[k]->alpha);
////			double alphfloor2 = floor(cell->cuts->vals[k2]->alpha);
////			double g01 = cell->cuts->vals[k]->alpha;
////			double g02 = cell->cuts->vals[k2]->alpha;
////			double gbar_0 = (cell->cuts->vals[k2]->alpha - cell->cuts->vals[k]->alpha)
////				- floor((cell->cuts->vals[k2]->alpha - cell->cuts->vals[k]->alpha));
////			if (gbar_0 > 0)
////			{
////				cut->beta[0] = 1;
////				cut->alpha = g01 + gbar_0 * (ceil(cell->cuts->vals[k2]->alpha - cell->cuts->vals[k]->alpha));
////
////				for (int v = 1; v < cell->master->mac; v++)
////				{
////
////					double betfloor1 = floor(cell->cuts->vals[k]->beta[v]);
////					double betfloor2 = floor(cell->cuts->vals[k2]->beta[v]);
////					double gj1 = cell->cuts->vals[k]->beta[v];
////					double gj2 = cell->cuts->vals[k2]->beta[v];
////					double gbar_j = (cell->cuts->vals[k2]->beta[v] - cell->cuts->vals[k]->beta[v]) -
////						floor(cell->cuts->vals[k2]->beta[v] - cell->cuts->vals[k]->beta[v]);
////					double beta1 = gbar_0 *(floor(cell->cuts->vals[k2]->beta[v] - cell->cuts->vals[k]->beta[v])) + gbar_j;
////					double beta2 = gbar_0 *(ceil(cell->cuts->vals[k2]->beta[v] - cell->cuts->vals[k]->beta[v]));
////
////
////					cut->beta[v] = min(beta1, beta2) + gj1;
////
////				}
////
////				bool repeatFlag = false;
////				if (cell->MIRcuts->vals == !NULL)
////				{
////					for (int pc = 0; pc < cell->MIRcuts->cnt; pc++)
////					{
////						if (cell->MIRcuts->vals[pc] != NULL && cut->alpha == cell->MIRcuts->vals[pc]->alpha &&
////							twoNorm(cut->beta, cell->MIRcuts->vals[pc]->beta, cell->master->mac) <= 0.001)
////						{
////							repeatFlag = true;
////							break;
////						}
////					}
////				}
////
////				if (repeatFlag == false)
////				{
////					cut->numSamples = cell->cuts->vals[k]->numSamples;
////					cell->MIRcuts->cnt++;
////					cutarr[cutnum] = cut;
////					cutnum++;
////				}
////
////			}
////
////		}
////	}
////#endif // defined(MIRSubbaddActive)
////
////
////	return cutarr;
////}//END formCut()
////
//// /*
////
//// form the pure integer fractional Chavatal-Gomory cuts
////
//// siavash tabrizian Apr 2020
////
//// */
////oneCut **purefracGMICut(probType **prob, cellType *cell, dVector Xvect, double lb) {
////	oneCut 	*cut;
////	oneCut 	**cutarr;
////	int 	k;
////
////	if (!(cutarr = (oneCut **) arr_alloc(cell->cuts->cnt + cell->cuts->cnt*cell->cuts->cnt / 2, oneCut *)))
////		errMsg("allocation", "formMIRcut", "cutarr", 0);
////
////	int cutnum = 0;
////	double alpha_coeff[1] = { 1.0 };
////
////	//// Creating frac GMI for benders cuts
////	for (int al = 0; al < 1; al++)
////	{
////		for (k = 0; k < cell->cuts->cnt; k++) {
////
////
////#if defined(GMIstreghtcuts1Active)
////			/* allocate memory to a new Chvatal Gomory cut */
////			cut = newCut(prob[0]->num->cols, 0, 0);
////
////			double bceil = ceil(alpha_coeff[al] * cell->cuts->vals[k]->alpha);
////			cut->beta[0] = bceil;
////
////			for (int v = 1; v < cell->master->mac; v++)
////			{
////
////				double aceil = ceil(alpha_coeff[al] * cell->cuts->vals[k]->beta[v]);
////				double maxpart = max(0, (bceil - aceil) / (1 + bceil));
////				cut->beta[v] = aceil - maxpart;
////
////			}
////
////			bool repeatFlag = false;
////			if (cell->GMIcuts->vals == !NULL)
////			{
////				for (int pc = 0; pc < cell->GMIcuts->cnt; pc++)
////				{
////					if (cell->GMIcuts->vals[pc] != NULL && cut->alpha == cell->GMIcuts->vals[pc]->alpha &&
////						twoNorm(cut->beta, cell->GMIcuts->vals[pc]->beta, cell->master->mac) <= 0.001)
////					{
////						repeatFlag = true;
////						break;
////					}
////				}
////			}
////
////			if (repeatFlag == false)
////			{
////				cut->numSamples = cell->cuts->vals[k]->numSamples;
////				cell->GMIcuts->cnt++;
////				cutarr[cutnum] = cut;
////				cutnum++;
////			}
////#else
////			/* allocate memory to a new Chvatal Gomory cut */
////			cut = newCut(prob[0]->num->cols, 0, 0);
////
////			double bceil = ceil(alpha_coeff[al] * cell->cuts->vals[k]->alpha);
////			cut->beta[0] = bceil;
////
////			for (int v = 1; v < cell->master->mac; v++)
////			{
////
////				double aceil = ceil(alpha_coeff[al] * cell->cuts->vals[k]->beta[v]);
////				cut->beta[v] = aceil;
////
////			}
////
////			bool repeatFlag = false;
////			if ( cell->GMIcuts->vals )
////			{
////				for (int pc = 0; pc < cell->GMIcuts->cnt; pc++)
////				{
////					if (cell->GMIcuts->vals[pc] != NULL && cut->alpha == cell->GMIcuts->vals[pc]->alpha &&
////						twoNorm(cut->beta, cell->GMIcuts->vals[pc]->beta, cell->master->mac) <= 0.001)
////					{
////						repeatFlag = true;
////						break;
////					}
////				}
////			}
////
////			if (repeatFlag == false)
////			{
////				cut->numSamples = cell->cuts->vals[k]->numSamples;
////				cell->GMIcuts->cnt++;
////				cutarr[cutnum] = cut;
////				cutnum++;
////			}
////#endif // defined(GMIstreghtcuts1Active)
////
////#if defined(GMIstreghtcuts2Active)
////			/* allocate memory to a new stregthen Chvatal Gomory cut */
////			/* based on Latchford and Lodi paper: Theorem 1*/
////			double fb = bceil - alpha_coeff[al] * cell->cuts->vals[k]->alpha;
////			cut = newCut(prob[0]->num->cols, 0, 0);
////
////			/* Find the max k */
////			int upk = ceil(1.0 / fb);
////			int maxk = 0;
////			for (int i = 1; i < upk; i++)
////			{
////				if ((1.0 / (float)(i + 1) <= fb) && (1.0 / (float)i > fb))
////				{
////					if (i > maxk)
////					{
////						maxk = i;
////					}
////				}
////			}
////			cut->beta[0] = (maxk + 1) * bceil;
////
////			/* Create N_p sets of indices*/
////			iVector N_p = (iVector)arr_alloc(maxk, double);
////			for (int v = 1; v < cell->master->mac; v++)
////			{
////
////				double aceil = ceil(alpha_coeff[al] * cell->cuts->vals[k]->beta[v]);
////				double fa = aceil - alpha_coeff[al] * cell->cuts->vals[k]->beta[v];
////				double ai = 0.0;
////				if (fa <= fb) ai += (maxk + 1) * aceil;
////				for (int p = 1; p <= maxk; p++)
////				{
////					if ((fb + ((p - 1)*(1 - fb)) / (float)maxk < fa) && (fa <= fb + (p*(1 - fb)) / (float)maxk))
////					{
////						ai += (maxk + 1)*aceil - p;
////					}
////				}
////				cut->beta[v] = ai;
////
////			}
////
////			repeatFlag = false;
////			if (cell->GMIcuts->vals == !NULL)
////			{
////				for (int pc = 0; pc < cell->GMIcuts->cnt; pc++)
////				{
////					if (cell->GMIcuts->vals[pc] != NULL && cut->alpha == cell->GMIcuts->vals[pc]->alpha &&
////						twoNorm(cut->beta, cell->GMIcuts->vals[pc]->beta, cell->master->mac) <= 0.001)
////					{
////						repeatFlag = true;
////						break;
////					}
////				}
////			}
////
////			if (repeatFlag == false)
////			{
////				cut->numSamples = cell->cuts->vals[k]->numSamples;
////				cell->GMIcuts->cnt++;
////				cutarr[cutnum] = cut;
////				cutnum++;
////			}
////#endif // defined(GMIstreghtcuts2Active)
////
////
////
////		}
////	}
////
////	return cutarr;
////}//END formCut()
//
//
