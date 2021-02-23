/*
 * smps.cpp
 *
 *  Created on: Mar 7, 2018
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "writer.hpp"

SMPSmodel::SMPSmodel() {

	/* Default parameters */
	numPeriods = 2;
	numStages = 2;

}//END constructor()


SMPSmodel::SMPSmodel(int t) {
    
    /* Default parameters */
    numPeriods = t;
    numStages = 2;
    
}//END constructor()

SMPSmodel::~SMPSmodel() {
 

}//END destructor()

