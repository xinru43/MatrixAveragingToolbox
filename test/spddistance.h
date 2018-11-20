/*
This is the test file for the Geometric mean of SPD problem defined in SPDMean.h and SPDMean.cpp.

---- WH
*/

#ifndef SPDDISTANCE_H
#define SPDDISTANCE_H

/*Help to debug the code*/
#include "ForDebug.h"

/*Output to console*/
#include <iostream>

/*Computational time*/
#include <ctime>
#include "Timer.h"

/*The global head file*/
#include "def.h"
#include "MyMatrix.h"


using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(SPDDISTANCE)
int main(void);
#endif

/*The main test function*/
void testspddistance(void);

#endif 
