/*
 * optimal.c
 *
 *  Created on: Apr 3, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "sdlp.h"

extern configType config;

BOOL optimal(probType **prob, cellType **cell, int T) {

	if ( cell[0]->k > config.MAX_ITER )
		return TRUE;

	return FALSE;
}//END optimal()
