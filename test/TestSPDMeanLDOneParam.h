/*
This is the test file for the Geometric mean of SPD problem defined in SPDMean.h and SPDMean.cpp.

---- WH
*/

#ifndef TESTSPDMEANSTEIN_H
#define TESTSPDMEANSTEIN_H

/*Help to debug the code*/
#include "ForDebug.h"

/*Output to console*/
#include <iostream>

/*Generate random number*/
#include "randgen.h"

/*Computational time*/
#include <ctime>

/*If this test file is called from Matlab, then functions in DriverMexProb.h are used.*/
#include "DriverMexProb.h"

/*Problem related classes*/
#include "Problem.h"
#include "SPDMeanLDOneParam.h"

/*Manifold related classes*/
#include "Manifold.h"
#include "SPDVector.h"
#include "SPDVariable.h"
#include "SPDManifold.h"

/*Linesearch based solvers*/
#include "SolversLS.h"
#include "RSD.h"
#include "RNewton.h"
#include "RCG.h"
#include "RBroydenFamily.h"
#include "RWRBFGS.h"
#include "RBFGS.h"
#include "LRBFGS.h"

/*Trust-region based solvers*/
#include "SolversTR.h"
#include "RTRSD.h"
#include "RTRNewton.h"
#include "RTRSR1.h"
#include "LRTRSR1.h"

/*The global head file*/
#include "def.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTSPDMEANLDONEPARAM)
int main(void);
#endif

/*The main test function*/
void testSPDMeanLDOneParam(void);

#endif // end of TESTSPDMEAN_H
