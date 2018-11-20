/*
This is the test file for the Geometric mean of SPD problem defined in SPDMean.h and SPDMean.cpp.

---- WH
*/

#ifndef TESTSPDEUCMEAN_H
#define TESTSPDEUCMEAN_H

/*Output to console*/
#include <iostream>

/*Generate random number*/
#include "Others/randgen.h"

/*Computational time*/
#include <ctime>

/*If this test file is called from Matlab, then functions in DriverMexProb.h are used.*/
#include "test/DriverMexProb.h"

/*Problem related classes*/
#include "Problems/Problem.h"
#include "Problems/SPDEucMean/SPDEucMean.h"

/*Manifold related classes*/
#include "Manifolds/Manifold.h"
#include "Manifolds/SPDEucManifold/SPDEucVector.h"
#include "Manifolds/SPDEucManifold/SPDEucVariable.h"
#include "Manifolds/SPDEucManifold/SPDEucManifold.h"

/*Linesearch based solvers*/
#include "Solvers/SolversLS.h"
#include "Solvers/RSD.h"
#include "Solvers/RNewton.h"
#include "Solvers/RCG.h"
#include "Solvers/RBroydenFamily.h"
#include "Solvers/RWRBFGS.h"
#include "Solvers/RBFGS.h"
#include "Solvers/LRBFGS.h"

/*Trust-region based solvers*/
#include "Solvers/SolversTR.h"
#include "Solvers/RTRSD.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/RTRSR1.h"
#include "Solvers/LRTRSR1.h"

/*The global head file*/
#include "Others/def.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTSPDEUCMEAN)
int main(void);
#endif

/*The main test function*/
void testSPDEucMean(void);

#endif // end of TESTSPDEucMEAN_H
