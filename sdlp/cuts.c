/*
 * cuts.c
 *
 *  Created on: Apr 2, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "sdlp.h"

void freeCutsType(cutsType *cuts) {
	int n;

	if (cuts->vals) {
		for ( n = 0; n < cuts->cnt; n++ )
			if ( cuts->vals[n]) freeOneCut(cuts->vals[n]);
		mem_free(cuts->vals);
	}
	mem_free(cuts);

}//END freeCutsType()

void freeOneCut(oneCut *cut) {

	if (cut->beta) mem_free(cut->beta);
	if (cut->iStar) mem_free(cut->iStar);
	mem_free(cut);

}//END freeOneCut()

