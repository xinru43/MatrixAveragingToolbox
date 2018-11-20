//
//
//  Test successfully
//  Cleaned the code compared with version 2.
//
//  Created by Xinru Yuan on 10/3/16.
//  Copyright 2016 Xinru Yuan. All rights reserved.
//

// #include "TestSPDRSD.h"
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

#include "TestSPDDLOneParamFP.h"

using namespace ROPTLIB;
// using namespace std;

double dis(double *X, integer n);
double distwo(double *X, double *Y, integer n);

#if !defined(MATLAB_MEX_FILE) && defined(TESTSPDDLONEPARAMFP)

int main()
{
     testSPDSPDDLOneParamFP();
    
#ifdef _WIN64
#ifdef _DEBUG
    _CrtDumpMemoryLeaks();
#endif
#endif
    return 0;
}

#endif





void testSPDSPDDLOneParamFP(void)
{
	printf("============= RSD: Using C++ =============\n");
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	tt = 0;
//	genrandseed(tt);

	/*Randomly generate a point on the SPD manifold*/
	integer n = 3, num = 3;
    double *initial_X = new double[n * n];
	for (integer i = 0; i < n; i++)
	{
		for (integer j = 0; j < n; j++)
		{
			initial_X[i + j * n] = 0;
		}
		initial_X[i + i * n] = 1;
	}

    
    integer N = n * n, length = n * n, Nnum = N * num, enditer, info;
    
    integer maxiter = 200;

	double *Rs = new double[n * n * num];
	double *Rstmp = new double[n * n];

    for (integer i = 0; i < num; i++)
	{
		for (integer j = 0; j < n * n; j++)
		Rstmp[j] = genrandnormal();// genrand_gaussian(); genrandnormal();
		dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, Rstmp, &n, Rstmp, &n, &GLOBAL::DZERO, Rs + i * n * n, &n);
	}
    
    double *AplusX = new double[N * num];
    double *iAplusX = new double[N * num];
    double *X2 = new double[n * n];
    double *tmpgrad = new double[n * n]; // tmpgrad is used to store gradient
    double *tmp = new double[N];
    double *tmp1 = new double[N];
    double *Ngrad = new double[N * num];
    
    
    //parameters
    double onehalf = 0.5, negone = -1.0;
    double ComTime;
    
    //make identity matrix
    double *I = new double[N];
    for(integer i = 0; i < N; i++) I[i] = 0.0;
    for(integer i = 0; i < n; i++) I[i*n+i] = 1.0;
    
    double *X1 = new double[N];
    dcopy_(&N, initial_X, &GLOBAL::IONE, X1, &GLOBAL::IONE);

    
    // Start the loop
    for (integer iter = 0; iter < maxiter; iter++)
    {
        for(integer h = 0; h < N; h++) X2[h] = 0.0;
        dcopy_(&Nnum, Rs, &GLOBAL::IONE, AplusX, &GLOBAL::IONE);
        
        for (integer i = 0; i < num; i++)
        {
            
            // Compute (Ai + X)/2 -> AplusX
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &onehalf, X1, &n, I, &n, &onehalf, AplusX + n * n * i, &n);
            
            // Cholesky factorization of AplusX -> AplusX
            dpotrf_(GLOBAL::L, &n, AplusX + n * n * i, &n, &info);
            for (integer j = 0; j < n; j++)
                for (integer k = j + 1; k < n; k++)
                    AplusX[j + k * n + i * n * n] = 0;
            
            // Compute inv(AplusX) by solving LY = I
            dcopy_(&length, I, &GLOBAL::IONE, iAplusX + n * n * i, &GLOBAL::IONE);
            dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &n, &n, AplusX + n * n * i, &n, iAplusX + n * n * i, &n, &info);
           
            // compute inv(AplusX) = inv(AplusX)^T * inv(AplusX)
            dcopy_(&length, iAplusX + n * n * i, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            dgemm_(GLOBAL::T, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, tmp, &n, tmp, &n, &GLOBAL::DZERO, iAplusX + n * n * i, &n);
           
            
            // X2 <- X2 + inv(AplusX)
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, iAplusX + n * n * i, &n, I, &n, &GLOBAL::DONE, X2, &n);
        }
        
        
        // compute inv(X2) -> tmp
        // Cholesky factorization of X2
        dpotrf_(GLOBAL::L, &n, X2, &n, &info);
        for (integer j = 0; j < n; j++)
            for (integer k = j + 1; k < n; k++)
                X2[j + k * n] = 0;
        
        // compute the inverse of chol(X2)
        dcopy_(&length, I, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &n, &n, X2, &n, tmp, &n, &info);
        

        // inv(tmp) = inv(tmp)^T * inv(tmp)
        dcopy_(&length, tmp, &GLOBAL::IONE, tmp1, &GLOBAL::IONE);
        dgemm_(GLOBAL::T, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, tmp1, &n, tmp1, &n, &GLOBAL::DZERO, tmp, &n);
               
        
        // compute num * tmp -> tmp
        for(integer i = 0; i < length; i++) tmp[i] = num * tmp[i];
        
        
        // ----------------- compute the norm of gradient g(tmp, tmp) -----------------
        for(integer i = 0; i < num; i++)
        {
            // compute X * iAplusX - I -> Ngrad
            dcopy_(&length, I, &GLOBAL::IONE, Ngrad + n * n * i, &GLOBAL::IONE);
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, X1, &n, iAplusX + n * n * i, &n, &negone, Ngrad + n * n * i, &n);
            
            // compute tmpgrad = sum(Ngrad)
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, Ngrad + n * n * i, &n, I, &n, &GLOBAL::DONE, tmpgrad, &n);
        }
        // compute gradient
        double gradscale = 1.0/(16 * num * num);
        dcopy_(&length, tmpgrad, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &gradscale, tmp, &n, tmp, &n, &GLOBAL::DONE, tmpgrad, &n);
        
//        gradSeries[iter] = 0.0;
//        for(integer j = 0; j < n; j++) gradSeries[iter] = gradSeries[iter] + tmpgrad[j*n+j];
//        gradSeries[iter] = std::sqrt(gradSeries[iter]);
//        
//        timeSeries[iter] = static_cast<double>(getTickCount() - starttime) / CLK_PS; // record the time for each iteration
//        
//        if(Debug >= 3 && iter > 0)
//        {
//            // compute distance between current iterate and identity matrix
//            disSeries[iter] = distwo(tmp, Xtrue, n);
//            truestepsizeSeries[iter] = distwo(X1, tmp, n);
//        }
        
        // update X1 <- tmp
        dcopy_(&length, tmp, &GLOBAL::IONE, X1, &GLOBAL::IONE);
        
    }
    
    delete[] Rstmp;
    delete[] Rs;
    delete[] initial_X;
    delete[] AplusX;
    delete[] iAplusX;
    delete[] X2;
    delete[] tmpgrad;
    delete[] Ngrad;
    delete[] tmp;
    delete[] tmp1;
    delete[] I;

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
    
    printf("============= RSD: Using mex =============\n");
    // deal with input
    double *RRs, *initial_X, *Xtrue;
    integer n, num, maxiter, Debug;
    double TOL, alpha;
    RRs = mxGetPr(prhs[0]);
    alpha = static_cast<double>(mxGetScalar(prhs[1]));
    initial_X = mxGetPr(prhs[2]);
    Xtrue = mxGetPr(prhs[3]);
    n = static_cast<integer>(mxGetScalar(prhs[4]));
    num = static_cast<integer>(mxGetScalar(prhs[5]));
    maxiter = static_cast<integer>(mxGetScalar(prhs[6]));
    TOL = static_cast<double>(mxGetScalar(prhs[7]));
    Debug = static_cast<integer>(mxGetScalar(prhs[8]));
    
    
    unsigned long starttime = getTickCount();
    CheckMemoryDeleted = new std::map<integer *, integer>;
    
    
    //parameters
    double onehalf = 0.5, negone = -1.0;
    double ComTime, err_g;
    double a1 = (1.0 - alpha)/2.0, a2 = (1.0 + alpha)/2.0;
    integer length = n * n, Nnum = length * num, enditer, info;
    
    //make identity matrix
    double *I = new double[length];
    for(integer i = 0; i < length; i++) I[i] = 0.0;
    for(integer i = 0; i < n; i++) I[i*n+i] = 1.0;
    
    double *X1 = new double[length];
    dcopy_(&length, initial_X, &GLOBAL::IONE, X1, &GLOBAL::IONE);
    
    
    double *timeSeries = new double[maxiter];
    double *gradSeries = new double[maxiter];
    double *disSeries = new double[maxiter];
    double *truestepsizeSeries = new double[maxiter];
    
     if(Debug >= 3)
     {
         // compute distance between current iterate and identity matrix
        disSeries[0] = distwo(X1, Xtrue, n);
     }

    
    double *AplusX = new double[Nnum];
    double *iAplusX = new double[Nnum];
    double *X2 = new double[length];
    double *tmpgrad = new double[length]; // tmpgrad is used to store gradient
    double *tmp = new double[length];
    double *tmp1 = new double[length];
    double *Ngrad = new double[Nnum];
  
    
// Start the loop
    for (integer iter = 0; iter < maxiter; iter++)
    {      
        for(integer h = 0; h < length; h++) X2[h] = 0.0;
        dcopy_(&Nnum, RRs, &GLOBAL::IONE, AplusX, &GLOBAL::IONE);
        
        for (integer i = 0; i < num; i++)
        {
          // Compute (Ai + X)/2 -> AplusX
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &a2, X1, &n, I, &n, &a1, AplusX + n * n * i, &n);

          // Cholesky factorization of AplusX -> AplusX
           dpotrf_(GLOBAL::L, &n, AplusX + n * n * i, &n, &info);
            for (integer j = 0; j < n; j++)
                for (integer k = j + 1; k < n; k++)
                    AplusX[j + k * n + i * n * n] = 0; 
           
          // Compute inv(AplusX) by solving LY = I
            dcopy_(&length, I, &GLOBAL::IONE, iAplusX + n * n * i, &GLOBAL::IONE);
            dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &n, &n, AplusX + n * n * i, &n, iAplusX + n * n * i, &n, &info);
            
            // compute inv(AplusX) = inv(AplusX)^T * inv(AplusX)
            dcopy_(&length, iAplusX + n * n * i, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            dgemm_(GLOBAL::T, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, tmp, &n, tmp, &n, &GLOBAL::DZERO, iAplusX + n * n * i, &n); 

            // X2 <- X2 + inv(AplusX)
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, iAplusX + n * n * i, &n, I, &n, &GLOBAL::DONE, X2, &n);
        }
        
        // compute inv(X2) -> tmp
        // Cholesky factorization of X2
        dpotrf_(GLOBAL::L, &n, X2, &n, &info);
        for (integer j = 0; j < n; j++)
            for (integer k = j + 1; k < n; k++)
                X2[j + k * n] = 0;
        
        // compute the inverse of chol(X2)
        dcopy_(&length, I, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &n, &n, X2, &n, tmp, &n, &info);
        

        // inv(tmp) = inv(tmp)^T * inv(tmp)
        dcopy_(&length, tmp, &GLOBAL::IONE, tmp1, &GLOBAL::IONE);
        dgemm_(GLOBAL::T, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, tmp1, &n, tmp1, &n, &GLOBAL::DZERO, tmp, &n);
                
        
        // compute num * tmp -> tmp  tmp is X2
        for(integer i = 0; i < length; i++) tmp[i] = num * tmp[i];

        
        // ----------------- compute the norm of gradient g(tmp, tmp) -----------------
        for(integer i = 0; i < length; i++) tmpgrad[i] = 0.0;
        for(integer i = 0; i < num; i++)
        {
            // compute X * iAplusX - I -> Ngrad
            dcopy_(&length, I, &GLOBAL::IONE, Ngrad + n * n * i, &GLOBAL::IONE);
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, X1, &n, iAplusX + n * n * i, &n, &negone, Ngrad + n * n * i, &n);
            
            // compute tmpgrad = sum(Ngrad)
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, Ngrad + n * n * i, &n, I, &n, &GLOBAL::DONE, tmpgrad, &n);
        } 
        // compute gradient
        dcopy_(&length, tmpgrad, &GLOBAL::IONE, tmp1, &GLOBAL::IONE);
        dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, tmp1, &n, tmp1, &n, &GLOBAL::DZERO, tmpgrad, &n);
        
        gradSeries[iter] = 0.0;
        for(integer j = 0; j < n; j++) gradSeries[iter] = gradSeries[iter] + tmpgrad[j*n+j];
        gradSeries[iter] = (2.0 * std::sqrt(gradSeries[iter]))/(1.0 - alpha);
        
        
        timeSeries[iter] = static_cast<double>(getTickCount() - starttime) / CLK_PS; // record the time for each iteration

        
        // compute distance between current iterate and identity matrix
          if(Debug >= 3 && iter > 0)
            {
                disSeries[iter] = distwo(tmp, Xtrue, n);
                truestepsizeSeries[iter] = distwo(X1, tmp, n);
            }
        
        
        // update X1 <- tmp
        dcopy_(&length, tmp, &GLOBAL::IONE, X1, &GLOBAL::IONE);


        // check stopping criteria
        err_g = gradSeries[iter]/gradSeries[0];
        
        if(err_g < TOL)
        {
            enditer = iter;
            break;
        }
        
    }
    
    ComTime = static_cast<double>(getTickCount() - starttime) / CLK_PS;

    // create output
    plhs[0] = mxCreateDoubleMatrix(maxiter, 1, mxREAL); // computation time
    plhs[1] = mxCreateDoubleScalar(ComTime); //  total computational time
    plhs[2] = mxCreateDoubleScalar(enditer); // number of iterations
    
    plhs[3] = mxCreateDoubleMatrix(length, 1, mxREAL); // X
    double *output3 = mxGetPr(plhs[3]);
    dcopy_(&length, X1, &GLOBAL::IONE, output3, &GLOBAL::IONE);
    
    plhs[4] = mxCreateDoubleMatrix(maxiter, 1, mxREAL); // distance
    plhs[5] = mxCreateDoubleMatrix(maxiter, 1, mxREAL); // true stepsize
    plhs[6] = mxCreateDoubleMatrix(maxiter, 1, mxREAL); // grad
    plhs[7] = mxCreateDoubleScalar(err_g); 
    
    double *plhstime = mxGetPr(plhs[0]);
    double *plhsdis = mxGetPr(plhs[4]);
    double *plhstruestepsize = mxGetPr(plhs[5]);
    double *plhsgrad = mxGetPr(plhs[6]);

    for(integer i = 0; i < maxiter; i++)
    {
        plhstime[i] = timeSeries[i];
        plhsdis[i] = disSeries[i];
        plhstruestepsize[i] = truestepsizeSeries[i];
        plhsgrad[i] = gradSeries[i];
//         printf("there: %.4e, %.4e, %.4e, %.4e\n", timeSeries[i], disSeries[i], stepsizeSeries[i], truestepsizeSeries[i]);
    }
      
    delete[] timeSeries;
    delete[] disSeries;  
    delete[] truestepsizeSeries;
    delete[] gradSeries;
    delete[] AplusX;
    delete[] iAplusX;
    delete[] X2;
    delete[] tmpgrad;
    delete[] Ngrad;
    delete[] tmp;
    delete[] I;
    delete[] tmp1;
    
    
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
    // compute distance between current iterate and the identity matrix
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







