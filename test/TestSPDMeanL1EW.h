/*
This is the test file for the Geometric mean of SPD problem defined in SPDMean.h and SPDMean.cpp.

---- WH
*/

#ifndef TESTSPDMEANL1EW_H
#define TESTSPDMEANL1EW_H

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
#include "SPDMean.h"
#include "SolversLS.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTSPDMEANL1EW)
int main(void);
#endif

/*The main test function*/
void testSPDMeanL1EW(void);

#endif 
