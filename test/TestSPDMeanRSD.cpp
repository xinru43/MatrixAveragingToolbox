//
//
//  Cleaned the code compared with version 2.
//
//  Created by Xinru Yuan on 10/3/16.
//  Copyright 2016 Xinru Yuan. All rights reserved.
//

// #include "TestSPDMeanRSD.h"
// 
// using namespace ROPTLIB;



/*Help to debug the code*/
// #include "ForDebug.h"
// #include <iostream>
// #include <ctime>
// #include "Timer.h"
// #include "def.h"
// #include "MyMatrix.h"
// #include "SPDMean.h"
// #include "SolversLS.h"

#include "TestSPDMeanRSD.h"

using namespace ROPTLIB;
// using namespace std;

double dis(double *X, integer n);
double distwo(double *X, double *Y, integer n);

#if !defined(MATLAB_MEX_FILE) && defined(TESTSPDMEANRSD)

int main()
{
     testSPDMeanRSD();
    
#ifdef _WIN64
#ifdef _DEBUG
    _CrtDumpMemoryLeaks();
#endif
#endif
    return 0;
}

#endif





void testSPDMeanRSD(void)
{
	printf("============= RSD: Using C++: needs to be implemented =============\n");
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	tt = 0;
//	genrandseed(tt);

	/*Randomly generate a point on the SPD manifold*/
	integer n = 100, num = 4;
    double *initial_X = new double[n * n];
	for (integer i = 0; i < n; i++)
	{
		for (integer j = 0; j < n; j++)
		{
			initial_X[i + j * n] = 0;
		}
		initial_X[i + i * n] = 1;
	}

    
    integer N = n * n;
    

	double *Rs = new double[n * n * num];
	double *Rstmp = new double[n * n];
	integer info;
	for (integer i = 0; i < num; i++)
	{
		for (integer j = 0; j < n * n; j++)
			Rstmp[j] = genrandnormal();
		dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, Rstmp, &n, Rstmp, &n, &GLOBAL::DZERO, Rs + i * n * n, &n);
	}
    delete[] Rstmp;
    delete[] Rs;
    delete[] initial_X;

};




#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

/*This function checks the number and formats of input parameters.
 * nlhs: the number of output in mxArray format
 * plhs: the output objects in mxArray format*
 * nrhs: the number of input in mxArray format
 * prhs: the input objects in mxArray format */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    printf("============= RSD =============\n");
    // input
    double *RRs, *initial_X, *Xtrue;
    integer n, num, maxiter, Retraction, StepType, Debug;
    double TOL;
    RRs = mxGetPr(prhs[0]);
    initial_X = mxGetPr(prhs[1]);
    Xtrue = mxGetPr(prhs[2]);
    n = static_cast<integer>(mxGetScalar(prhs[3]));
    num = static_cast<integer>(mxGetScalar(prhs[4]));
    maxiter = static_cast<integer>(mxGetScalar(prhs[5]));
    TOL = static_cast<double>(mxGetScalar(prhs[6]));
    Retraction = static_cast<integer>(mxGetScalar(prhs[7]));
    StepType = static_cast<integer>(mxGetScalar(prhs[8]));
    Debug = static_cast<integer>(mxGetScalar(prhs[9]));
    
    
    unsigned long starttime = getTickCount();
    CheckMemoryDeleted = new std::map<integer *, integer>;
    
    
    //parameters
    double a, stepsize, ch, Delta, onehalf = 0.5;
    double ComTime;
    integer N = n * n, enditer, info;
    
    
    double *X1 = new double[N];
    dcopy_(&N, initial_X, &GLOBAL::IONE, X1, &GLOBAL::IONE);
    
    double *Rs = new double[N * num];
    integer Nnum = N * num;
    dcopy_(&Nnum, RRs, &GLOBAL::IONE, Rs, &GLOBAL::IONE);
        
    
    double *timeSeries = new double[maxiter];
    double *gradSeries = new double[maxiter];
    double *disSeries = new double[maxiter];
    double *truestepsizeSeries = new double[maxiter];
    double *stepsizeSeries = new double[maxiter];


    double *Lx = new double[N]; // store the Cholesky factorization of each iterate X
    double *Z = new double[N * num];  // store R0^{-1}R
    double *S = new double[N];  // store the gradient
    double *LiStmp = new double[N];
    double *Y = new double[N];
    double *tmp = new double[N];
    double *eigenvalues = new double[n + 2 * N];
    double *eigenvectors = eigenvalues + n;
    double *eigenvectorsD = eigenvectors + N;
    Matrix E(eigenvalues, n, 1), V(eigenvectors, n, n), VD(eigenvectorsD, n, n);

    
    //make identity matrix
    double *I = new double[N];
    for(integer i = 0; i < N; i++) I[i] = 0.0;
    for(integer i = 0; i < n; i++) I[i*n+i] = 1.0;

    
    integer iter = 0;
    
    
     if(Debug >= 3)
     {
         // compute distance between current iterate and the true solution
        disSeries[iter] = distwo(X1, Xtrue, n);
     }
    
         
// Start the loop
    for (integer iter = 0; iter < maxiter; iter++)
    {
        
        // get cholesky of X0 and stored in Lx
        dcopy_(&N, X1, &GLOBAL::IONE, Lx, &GLOBAL::IONE);
        dpotrf_(GLOBAL::L, &n, Lx, &n, &info);
        for (integer j = 0; j < n; j++)
        {
            for (integer k = j + 1; k < n; k++)
            {
                Lx[j + k * n] = 0; //manually make the upper triangle 0: X0 = X0 * X0^T
            }
        }
        

        Delta = 0.0;
        for(integer h = 0; h < N; h++) S[h] = 0.0; // S is used to store the gradient
        
        
        for (integer i = 0; i < num; i++)
        {
            // solve linear system R0*Z=Rs, then Z=R0^{-1}Rs
            dcopy_(&N, Rs + i * N, &GLOBAL::IONE, Z + i * N, &GLOBAL::IONE);
            dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &n, &n, Lx, &n, Z + i * N, &n, &info);
            
            // compute Z{i} * Z{i}' and store in Z{i}
            dcopy_(&N, Z + i * N, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, tmp, &n, tmp, &n, &GLOBAL::DZERO, Z + i * N, &n);
            
            // eigenvalue decomposition Z{i}
            Matrix MMt(Z + i * N, n, n);
            Matrix::EigenSymmetricM(GLOBAL::L, MMt, E, V);
            
            // compute the condition number of Z
            ch = eigenvalues[n-1]/eigenvalues[0];
            
            Delta = Delta + ((ch + 1.0) * log(ch))/(ch - 1.0);
            
            //compute Z <- log(Z'*Z)
            dcopy_(&N, eigenvectors, &GLOBAL::IONE, VD.matrix, &GLOBAL::IONE);
            for (integer i = 0; i < n; i++)
            {
                a = log(eigenvalues[i]);
                dscal_(&n, &a, eigenvectors + i * n, &GLOBAL::IONE);
            }
            Matrix MMt2(Z + i * N, n, n);
            Matrix::DGEMM(GLOBAL::DONE, V, false, VD, true, GLOBAL::DZERO, MMt2);
            
            // Z <- (Z + Z')/2
            dcopy_(&N, Z + i * N, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            dgemm_(GLOBAL::T, GLOBAL::N, &n, &n, &n, &onehalf, tmp, &n, I, &n, &onehalf, Z + i * N, &n);
            
            // S = S + Z
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, Z + i * N, &n, I, &n, &GLOBAL::DONE, S, &n);
        }        
        
        
        // compute stepsize
        if (StepType == 1) 
        {
            stepsize = 2.0/Delta;  // RL
        }
        else if (StepType == 0)
        {
            stepsize = 4.0/(2.0 * num + Delta); //QR
        }
        
        stepsizeSeries[iter] = stepsize;
        
        
        // compute Xk <- R0 * S * R0'//
        // S <- R0 * S
        dcopy_(&N, S, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, Lx, &n, tmp, &n, &GLOBAL::DZERO, S, &n);
        
        // S <- S * R0'
        dcopy_(&N, S, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, tmp, &n, Lx, &n, &GLOBAL::DZERO, S, &n);
        
        
        // ----------------- compute the norm of gradient
        // compute |S| = tr(S * X^{-1} * S * X^{-1})
        // compute Y1 <- L^{-1}S by solving L*Y1 = S
        dcopy_(&N, S, &GLOBAL::IONE, Y, &GLOBAL::IONE);
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &n, &n, Lx, &n, Y, &n, &info);
        
        // transpose Y1
        dcopy_(&N, Y, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::T, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, tmp, &n, I, &n, &GLOBAL::DZERO, Y, &n);
        
        // compute L * Y = Y1
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &n, &n, Lx, &n, Y, &n, &info);
        
        // compute Y <- Y * Y
        dcopy_(&N, Y, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, tmp, &n, tmp, &n, &GLOBAL::DZERO, Y, &n);
        
        gradSeries[iter] = 0.0;
        for(integer j = 0; j < n; j++) gradSeries[iter] = gradSeries[iter] + Y[j*n+j];
        gradSeries[iter] = std::sqrt(gradSeries[iter]);
        
//        printf("iter: %d, |grad|: %.4e\n", iter, gradSeries[iter]);
        
        timeSeries[iter] = static_cast<double>(getTickCount() - starttime) / CLK_PS; // record the time for each iteration                
        
        
        // S <- theta * S
        dcopy_(&N, S, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &stepsize, tmp, &n, I, &n, &GLOBAL::DZERO, S, &n);
        

        //Retraction:  S <- X0 + S
        if (Retraction == 1)
        {
          dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, X1, &n, I, &n, &GLOBAL::DONE, S, &n);
        }
        else if(Retraction == 2)
        {
            // X2 <- X1 + eta + 0.5 * eta * X^{-1} * eta
            dcopy_(&N, S, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            // Solve the linear system L X = E, i.e., X = Lx^{-1} eta. The solution X is stored in LiE
            dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &n, &n, Lx, &n, tmp, &n, &info);
            
            //Compute X2 = LiE^T LiE = E L^{-T} L^{-1} E
            dgemm_(GLOBAL::T, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, tmp, &n, tmp, &n, &GLOBAL::DZERO, LiStmp, &n);
            
            /*X2 <-- 0.5 * X2*/
            dscal_(&N, &onehalf, LiStmp, &GLOBAL::IONE);
            
            /*S <-- S + LiStmp */
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, LiStmp, &n, I, &n, &GLOBAL::DONE, S, &n);
            
            /*S <-- x + S*/
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, X1, &n, I, &n, &GLOBAL::DONE, S, &n);
        }
                
            if(Debug >= 3 && iter > 0)
            {
                // compute distance between current iterate and identity matrix
                disSeries[iter] = distwo(X1, Xtrue, n);
                truestepsizeSeries[iter] = distwo(X1, S, n);
            }
        
//        printf("dis: %.4e\n", disSeries[iter]);
  
        
        // Update X_old <- X_new: S is the new iterate, initial_X is the old iterate
        dcopy_(&N, S, &GLOBAL::IONE, X1, &GLOBAL::IONE);

        // check stopping criteria
        if(gradSeries[iter] < TOL)
        {
            enditer = iter;
            break;
        }
        
    }
    
    ComTime = static_cast<double>(getTickCount() - starttime) / CLK_PS;

    

    // create output
    plhs[0] = mxCreateDoubleMatrix(maxiter, 1, mxREAL); // computation time
    plhs[1] = mxCreateDoubleScalar(ComTime); //  total computational time
    plhs[2] = mxCreateDoubleScalar(maxiter); // number of iterations
    
    plhs[3] = mxCreateDoubleMatrix(N, 1, mxREAL); // X
    double *output3 = mxGetPr(plhs[3]);
    dcopy_(&N, initial_X, &GLOBAL::IONE, output3, &GLOBAL::IONE);
    
    plhs[4] = mxCreateDoubleMatrix(maxiter, 1, mxREAL); // distance
    plhs[5] = mxCreateDoubleMatrix(maxiter, 1, mxREAL); // stepsize
    plhs[6] = mxCreateDoubleMatrix(maxiter, 1, mxREAL); // true stepsize
    plhs[7] = mxCreateDoubleMatrix(maxiter, 1, mxREAL); // grad
    
    
    double *plhstime = mxGetPr(plhs[0]);
    double *plhsdis = mxGetPr(plhs[4]);
    double *plhsstepsize = mxGetPr(plhs[5]);
    double *plhstruestepsize = mxGetPr(plhs[6]);
    double *plhsgrad = mxGetPr(plhs[7]);

    for(integer i = 0; i < maxiter; i++)
    {
        plhstime[i] = timeSeries[i];
        plhsdis[i] = disSeries[i];
        plhsstepsize[i] = stepsizeSeries[i];
        plhstruestepsize[i] = truestepsizeSeries[i];
        plhsgrad[i] = gradSeries[i];
//         printf("there: %.4e, %.4e, %.4e, %.4e\n", timeSeries[i], disSeries[i], stepsizeSeries[i], truestepsizeSeries[i]);
    }
      
    delete[] timeSeries;
    delete[] disSeries;  
    delete[] stepsizeSeries;
    delete[] truestepsizeSeries;
    delete[] gradSeries;
    delete[] Z;
    delete[] S;
    delete[] X1;
    delete[] Rs;
    delete[] I;
    delete[] LiStmp;
    delete[] Lx;
    delete[] Y;
    
    
    std::map<integer *, integer>::iterator iteriter = CheckMemoryDeleted->begin();
	for (iteriter = CheckMemoryDeleted->begin(); iteriter != CheckMemoryDeleted->end(); iteriter++)
	{
		if (iteriter->second != 1)
            printf("Global address: %d", iteriter->first);
            printf(", shared times: %d\n", iteriter->second);
	}
	delete CheckMemoryDeleted;
    return;
    
}

#endif




double dis(double *X, integer n)
{
    // compute distance between current iterate and identity matrix
    double *eigenvaluesXnew = new double[n + n * n];
    double *eigenvectorsXnew = eigenvaluesXnew + n;
    double *tmp = new double[n * n];
    
    double disSeries;
    integer N = n * n;
    
    dcopy_(&N, X, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
    Matrix Enew(eigenvaluesXnew, n, 1), Vnew(eigenvectorsXnew, n, n);
    Matrix MMtnew(tmp, n, n);
    Matrix::EigenSymmetricM(GLOBAL::L, MMtnew, Enew, Vnew);
    disSeries = 0.0;
    for (integer ii = 0; ii < n; ii++)
    {
        disSeries = disSeries+log(eigenvaluesXnew[ii])*log(eigenvaluesXnew[ii]);
    }
    disSeries = sqrt(disSeries);
    delete[] eigenvaluesXnew;
    delete[] tmp;
    return disSeries;
    
};




double distwo(double *X, double *Y, integer n)
{
    double *Xold = new double[n * n];
    double *Xnew = new double[n * n];
    double *tmp = new double[n * n];
    double *eigenvaluesXold = new double[n + n * n];
    double *eigenvectorsXold = eigenvaluesXold + n;
    double *eigenvaluesXnew = new double[n + n * n];
    double *eigenvectorsXnew = eigenvaluesXnew + n;
    integer N = n * n;
    integer info;
    double truestepsize;
    
    dcopy_(&N, X, &GLOBAL::IONE, Xold, &GLOBAL::IONE);
    dcopy_(&N, Y, &GLOBAL::IONE, Xnew, &GLOBAL::IONE);
    
    // get cholesky of Xnew and Xold
    dpotrf_(GLOBAL::L, &n, Xnew, &n, &info);
    dpotrf_(GLOBAL::L, &n, Xold, &n, &info);
    for (integer j = 0; j < n; j++)
    {
        for (integer k = j + 1; k < n; k++)
        {
            Xnew[j + k * n] = 0; //manually make the upper triangle 0: X0 = X0 * X0^T
            Xold[j + k * n] = 0;
        }
    }
    // Solve linear system Xnew * W = Xold, i.e. W = Xnew^{-1}*Xold -> Xold
    dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &n, &n, Xnew, &n, Xold, &n, &info);
    
    // Compute Xold <- Xold * Xold'
    dcopy_(&N, Xold, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
    dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, tmp, &n, tmp, &n, &GLOBAL::DZERO, Xold, &n);
    
    // eigenvalue decomposition of Xold
    Matrix Eold(eigenvaluesXold, n, 1), Vold(eigenvectorsXold, n, n);
    Matrix MMtold(Xold, n, n);
    Matrix::EigenSymmetricM(GLOBAL::L, MMtold, Eold, Vold);
    
    truestepsize = 0.0;
    for (integer ii = 0; ii < n; ii++)
    {
        truestepsize = truestepsize+log(eigenvaluesXold[ii])*log(eigenvaluesXold[ii]);
    }
    truestepsize = sqrt(truestepsize);
    
    delete[] Xold;
    delete[] Xnew;
    delete[] tmp;
    delete[] eigenvaluesXold;
    delete[] eigenvaluesXnew;
    
    return truestepsize;
};







