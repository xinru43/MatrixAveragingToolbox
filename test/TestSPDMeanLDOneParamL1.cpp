
#include "TestSPDMeanLDOneParamL1.h"

using namespace ROPTLIB;

/*If the file is not compiled in Matlab and TESTSTIEBROCKETT is defined in def.h file, then using the following
main() function as the entrance. */
#if !defined(MATLAB_MEX_FILE) && defined(TESTSPDMEANLDONEPARAML1)

int main(void)
{
	testSPDMeanLDOneParamL1();

#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}

#endif

void testSPDMeanLDOneParamL1(void)
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	tt = 0;
	genrandseed(tt);

	/*Randomly generate a point on the SPD manifold*/
	integer n = 100, num = 4;
    double alpha = 0.5;
	SPDVariable SPDX(n);
	double *initialX = SPDX.ObtainWriteEntireData();
	for (integer i = 0; i < n; i++)
	{
		for (integer j = 0; j < n; j++)
		{
			initialX[i + j * n] = 0;
		}
		initialX[i + i * n] = 1;
	}

	// Define the manifold
	SPDManifold Domain(n);
	Domain.SetHasHHR(true); /*set whether the manifold uses the idea in [HGA2015, Section 4.3] or not*/

	double *Ls = new double[n * n * num + n * n];
    double *LDAs = new double[num];
    double *Xtrue = new double[n * n];
	double *tmp = Ls + n * n * num;
	integer info;
	for (integer i = 0; i < num; i++)
	{
		for (integer j = 0; j < n * n; j++)
			tmp[j] = genrandnormal();

		dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, tmp, &n, tmp, &n, &GLOBAL::DZERO, Ls + i * n * n, &n);
    }
    
    for(integer i = 0; i < num; i++) LDAs[i] = 0.0; // <--------------- needs to modify
    
    
    
	// Define the problem
	SPDMeanLDOneParamL1 Prob(Ls, LDAs, n, num, alpha);
	/*The domain of the problem is a SPD manifold*/
	Prob.SetDomain(&Domain);

	//Prob.CheckGradHessian(&SPDX);

	/*Output the parameters of the domain manifold*/
	Domain.CheckParams();

	/*Check the correctness of the manifold operations*/
	//Domain.CheckIntrExtr(&SPDX);
	//Domain.CheckRetraction(&SPDX);
	//Domain.CheckDiffRetraction(&SPDX);
	//Domain.CheckLockingCondition(&SPDX);
	//Domain.CheckcoTangentVector(&SPDX);
	//Domain.CheckIsometryofVectorTransport(&SPDX);
	//Domain.CheckIsometryofInvVectorTransport(&SPDX);
	//Domain.CheckVecTranComposeInverseVecTran(&SPDX);
	//Domain.CheckTranHInvTran(&SPDX);
	//Domain.CheckHaddScaledRank1OPE(&SPDX);

	// test LRBFGS
	std::cout << "********************************Test Geometric mean in LRBFGS*************************************" << std::endl;
	LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &SPDX);
	LRBFGSsolver->LineSearch_LS = ARMIJO;
	LRBFGSsolver->Debug = ITERRESULT; //ITERRESULT;// 
	LRBFGSsolver->Max_Iteration = 20;
	LRBFGSsolver->Tolerance = 1e-10;
	LRBFGSsolver->Accuracy = 1e-4;
	LRBFGSsolver->Finalstepsize = 1;
	LRBFGSsolver->CheckParams();
	LRBFGSsolver->Run();
	delete LRBFGSsolver;

	delete[] Ls;
};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 4)
	{
		mexErrMsgTxt("The number of arguments should be at least four.\n");
	}
	double *Ls, *LDAs, *X, *Xopt, *soln;
    double alpha;
    
    Ls = mxGetPr(prhs[0]);
    LDAs = mxGetPr(prhs[1]);
    alpha = static_cast<double> (mxGetScalar(prhs[2]));
    X = mxGetPr(prhs[3]);
    
    /* dimensions of input matrices */ 
    integer n, N, HasHHR;
    n = mxGetM(prhs[3]);
    
	const mwSize *ptrdims = mxGetDimensions(prhs[0]);
	if (mxGetNumberOfDimensions(prhs[0]) == 2)
		N = 1;
	else
		N = ptrdims[2];    

	if (ptrdims[1] != n || ptrdims[0] != n)
	{
		mexErrMsgTxt("The size of matrix C is not correct.\n");
	}
	if (mxGetM(prhs[3]) != n || mxGetN(prhs[3]) != n)
	{
		mexErrMsgTxt("The size of the initial X is not correct!\n");
	}
    
    HasHHR = static_cast<integer> (mxGetScalar(prhs[4]));

    // create output matrix
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
	Xopt = mxGetPr(plhs[0]);

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
    
    
	// Obtain an initial iterate from the input
	SPDVariable *SPDX = new SPDVariable(n);
	double *initialX = SPDX->ObtainWriteEntireData();
	for (integer i = 0; i < n * n; i++)
	{
		initialX[i] = X[i];
	}
    
    
    SPDVariable *SPDsoln = nullptr;
    if (nrhs >= 7)
	{
		double *soln = mxGetPr(prhs[6]);
		SPDsoln = new SPDVariable(n);
		double *SPDsolnptr = SPDsoln->ObtainWriteEntireData();
		for (integer i = 0; i < n * n; i++)
		{
			SPDsolnptr[i] = soln[i];
		}
	}


	// Define the manifold
	SPDManifold Domain(n);
	// Define the SPDMean problem
	SPDMeanLDOneParamL1 Prob(Ls, LDAs, n, N, alpha);
	Prob.SetDomain(&Domain);
	Domain.SetHasHHR((HasHHR != 0));
	ParseSolverParamsAndOptimizing(prhs[5], &Prob, SPDX, SPDsoln, plhs);
    
	//Domain.CheckParams();
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			std::cout << "Global address:" << iter->first << ", sharedtimes:" << iter->second << std::endl;
	}
	delete CheckMemoryDeleted;
	delete SPDX;
    delete SPDsoln;
	return;
}

#endif
