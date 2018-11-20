//  Created by Xinru Yuan on 10/3/16.
//  Copyright 2016 Xinru Yuan. All rights reserved.

/*Help to debug the code*/
// #include "ForDebug.h"
// #include <iostream>
// #include <ctime>
// #include "Timer.h"
// #include "def.h"
// #include "MyMatrix.h"
// #include "SPDMean.h"
// #include "SolversLS.h"

#include "TestSPDMeanLDOneParamL1EW.h"

using namespace ROPTLIB;

double dis(double *X, integer n);
double distwo(double *X, double *Y, integer n);

#if !defined(MATLAB_MEX_FILE) && defined(TESTSPDMEANLDONEPARAML1EW)

int main()
{
    testSPDMeanL1LDOneParamEW();
    
#ifdef _WIN64
#ifdef _DEBUG
    _CrtDumpMemoryLeaks();
#endif
#endif
    return 0;
}

#endif


void testSPDMeanLDOneParamL1EW(void)
{
    printf("============= LogDet alpha-divergence median computation by EW method =============\n");
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
    
    
    integer N = n, length = n * n, enditer, info;
    integer maxiter = 500;
    
    
    double *Ls = new double[n * n * num];
    double *Rstmp = new double[n * n];
    
    for (integer i = 0; i < num; i++)
    {
        for (integer j = 0; j < n * n; j++)
            Rstmp[j] = genrandnormal();// genrand_gaussian(); genrandnormal();
        dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, Rstmp, &n, Rstmp, &n, &GLOBAL::DZERO, Ls + i * n * n, &n);
        
        
        dpotrf_(GLOBAL::L, &n, Ls + n * n * i, &n, &info);
        for (integer j = 0; j < n; j++)
            for (integer k = j + 1; k < n; k++)
                Ls[j + k * n + i * n * n] = 0;
        
    }
    
    
    unsigned long starttime = getTickCount();
    
    double minval, negone = -1.0, t;
    double ComTime, err_g;
    integer idx;
    
    double *X1 = new double[n * n];
    dcopy_(&length, initial_X, &GLOBAL::IONE, X1, &GLOBAL::IONE);
    
    double *timeSeries = new double[maxiter];
    double *funsSeries = new double[maxiter];
    double *gradSeries = new double[maxiter];
    double *disSeries = new double[maxiter];
    
    double *logLXL = new double[n * n * num];
    double *coef = new double[num];
    double *eta = new double[n * n];
    double *Lx = new double[n * n];
    double *tmp = new double[n * n];
    double *X2 = new double[n * n];
    
    
// Start the loop
    for (integer iter = 0; iter < maxiter; iter++)
    {
//         Cholesky factorization of X1, and store in Lx
        dcopy_(&length, X1, &GLOBAL::IONE, Lx, &GLOBAL::IONE);
        
        dpotrf_(GLOBAL::L, &n, Lx, &n, &info);
        for (integer j = 0; j < n; j++)
            for (integer k = j + 1; k < n; k++)
                Lx[j + k * n] = 0;
        
//   --------------  distance from A_i to X_k  --------------
        for (integer i = 0; i < num; i++)
        {
            dcopy_(&length, Lx, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            /*Solve the linear system Li X = Lx, i.e., X = Li^{-1} Lx. The solution X is stored in LiiLx.
             * Note that Li is a lower triangular matrix.
             * Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
            dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, Ls + n * n * i, &N, tmp, &N, &info);
            if (info != 0)
            {
                printf("The cholesky decompsotion in SPDMean::f failed with info:%d!\n", info);
            }
            dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, tmp, &N, tmp, &N, &GLOBAL::DZERO, logLXL + n * n * i, &N);
            Matrix MMt(logLXL + n * n * i, n, n);
            Matrix::LogSymmetricM(GLOBAL::L, MMt, MMt);
            coef[i] = dnrm2_(&length, logLXL + n * n * i, &GLOBAL::IONE);
        }
        
        
        minval = coef[0];
        idx = 0;
        for (integer i = 1; i < num; i++)
        {
            if(coef[i] < minval)
            {
                minval = coef[i];
                idx = i;
            }
        }
        
        
        
        funsSeries[iter] = 0;
        for (integer i = 0; i < num; i++) funsSeries[iter] = funsSeries[iter] + coef[i];
        
//          --------------  compute eta  --------------
//          logLXL * L^T -> tmp
        dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, logLXL + n * n * idx, &N, Ls + n * n * idx, &N, &GLOBAL::DZERO, tmp, &N);
//         ForDebug::Print("tmp", tmp, n, n);
        
//          L^{-T} * tmp -> tmp, i.e.,  L^{-T} * logLXL * L^T -> tmp
        dtrtrs_(GLOBAL::L, GLOBAL::T, GLOBAL::N, &N, &N, Ls + n * n * idx, &N, tmp, &N, &info);
        
        if (info != 0)
        {
            printf("The cholesky decompsotion in SPDMean::RieGrad failed with info:%d!\n", info);
        }
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, X1, &N, tmp, &N, &GLOBAL::DZERO, eta, &N);
        
        
//      eta -> -eta
        dscal_(&length, &negone, eta, &GLOBAL::IONE);
        
//          --------------  update next iterate  --------------
        t = 1.0/(1.0 + iter);
        dscal_(&length, &t, eta, &GLOBAL::IONE);
        
//      Cholesky factorization of eta, and stored in eta
        dpotrf_(GLOBAL::L, &n, eta, &n, &info);
        for (integer j = 0; j < n; j++)
            for (integer k = j + 1; k < n; k++)
                eta[j + k * n] = 0;
        
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, Ls + n * n * idx, &N, eta, &N, &info);
        if (info != 0)
        {
            printf("The cholesky decompsotion in SPDMean::f failed with info:%d!\n", info);
        }
        dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, eta, &N, eta, &N, &GLOBAL::DZERO, X2, &N);
        Matrix MMt(X2, n, n);
        Matrix::ExpSymmetricM(GLOBAL::L, MMt, MMt);
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, Lx, &N, X2, &N, &GLOBAL::DONE, X2, &N);
        dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, X2, &N, Lx, &N, &GLOBAL::DONE, X2, &N);
        
        
        //         update X1 <- X2
        dcopy_(&length, X2, &GLOBAL::IONE, X1, &GLOBAL::IONE);
//         ForDebug::Print("X2", X2, n, n);
        
        
        
        timeSeries[iter] = static_cast<double>(getTickCount() - starttime) / CLK_PS; // record the time for each iteration
        
        // check stopping criteria
    }
    
    
    ComTime = static_cast<double>(getTickCount() - starttime) / CLK_PS;
    
    
    delete[] timeSeries;
    delete[] disSeries;
    delete[] gradSeries;
    delete[] funsSeries;
    
    delete[] logLXL;
    delete[] coef;
    delete[] eta;
    delete[] Lx;
    delete[] tmp;
    delete[] X1;
    delete[] X2;
    delete[] Rstmp;
    
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
    
    printf("============= LogDet alpha-divergence median computation by EW method =============\n");
    // deal with input
    double *Ls, *LDAs, *initial_X, *Xtrue;
    integer n, num, maxiter, Debug;
    double TOL, alpha, stepsize;
    Ls = mxGetPr(prhs[0]);
    LDAs = mxGetPr(prhs[1]);
    alpha = static_cast<double>(mxGetScalar(prhs[2]));
    initial_X = mxGetPr(prhs[3]);
    Xtrue = mxGetPr(prhs[4]);
    n = static_cast<integer>(mxGetScalar(prhs[5]));
    num = static_cast<integer>(mxGetScalar(prhs[6]));
    maxiter = static_cast<integer>(mxGetScalar(prhs[7]));
    stepsize = static_cast<double>(mxGetScalar(prhs[8]));
    TOL = static_cast<double>(mxGetScalar(prhs[9]));
    Debug = static_cast<integer>(mxGetScalar(prhs[10]));
    
    
//     time start
    unsigned long starttime = getTickCount();
    CheckMemoryDeleted = new std::map<integer *, integer>;
    
    
    integer length = n * n, N = n, Nnum = n * n * num, info, enditer;
    double minval, negone = -1.0, t, coefa, scale;
    double ComTime, err_g, sigma;
    double a1 = (1.0 - alpha)/2.0, a2 = (1.0 + alpha)/2.0;
    integer idx;
    
    // X1 is used to store initial_X
    double *X1 = new double[n * n];
    dcopy_(&length, initial_X, &GLOBAL::IONE, X1, &GLOBAL::IONE);
    
    double *timeSeries = new double[maxiter];
    double *funsSeries = new double[maxiter];
    double *gradSeries = new double[maxiter];
    double *disSeries = new double[maxiter];
    
    double *coef = new double[num];
    double *eta = new double[n * n];
    double *L = new double[n * n];
    double *tmp = new double[n * n];
    double *X2 = new double[n * n];
    
    double *CLplusX = new double[n * n * num];
    double *LLplusX = new double[n * n * num];
    double *iLplusX = new double[n * n * num];
    
    double *I = new double[n * n]; // make identity matrix
    for (integer i = 0; i < length; i++) I[i] = 0.0;
    for (integer i = 0; i < n; i++) I[i*n+i] = 1.0;
    
    
    if(Debug >= 3)
    {
        // compute distance between current iterate and identity matrix
        disSeries[0] = distwo(X1, Xtrue, n);
        gradSeries[0] = 0.0;
    }
    
// Start the loop
    for (integer iter = 0; iter < maxiter; iter++)
    {
        
        dcopy_(&Nnum, Ls, &GLOBAL::IONE, CLplusX, &GLOBAL::IONE);
        
//   --------------  divergence from A_i to X_k  --------------
        for (integer i = 0; i < num; i++)
        {
            // compute a1 * Ai + a2 * x
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &a2, X1, &N, I, &N, &a1, CLplusX + n * n * i, &N);
            
            // Cholesky factorization of CLplusX
            dpotrf_(GLOBAL::L, &N, CLplusX + n * n * i, &N, &info);
            for (integer j = 0; j < n; j++)
                for (integer k = j + 1; k < n; k++)
                    CLplusX[j + k * n + i * n * n] = 0;
            
            // compute the determinant of CLplusX and stored in coef
            coef[i] = 1.0;
            for (integer j = 0; j < n; j++)
            {
                coef[i] = coef[i] * CLplusX[j * n + j + i * n * n];
            }
            coef[i] = log(coef[i] * coef[i]);
        }
        
        // compute the determinant of X: product of diagonal elements, stored in L
        dcopy_(&length, X1, &GLOBAL::IONE, L, &GLOBAL::IONE);
        
        dpotrf_(GLOBAL::L, &n, L, &n, &info);
        for (integer j = 0; j < n; j++)
            for (integer k = j + 1; k < n; k++)
                L[j + k * n] = 0;
        
        double detX = 1.0;
        for (integer j = 0; j < n; j++)
        {
            detX = detX * L[j * n + j];
        }
        detX = a2 * log(detX * detX);
        
        
        for (integer i = 0; i < num; i++)
        {
            coef[i] = coef[i] - detX - a1 * LDAs[i];
        }
        
        
//   --------------  check if div(Xk, Ai) == 0: find the smallest div(Xk, Ai)  --------------
        minval = coef[0];
        idx = 0;
        for (integer i = 1; i < num; i++)
        {
            if(coef[i] < minval)
            {
                minval = coef[i];
                idx = i;
            }
        }
        
        
//          --------------  compute eta  --------------
        dcopy_(&Nnum, CLplusX, &GLOBAL::IONE, LLplusX, &GLOBAL::IONE);
        
        if (minval > 1e-16)
        {
            // compute sigma = sum 1/dist(Xk, Ai)
            sigma = 0.0;
            for (integer i = 0; i < num; i++)
            {
                sigma = sigma + 1.0/coef[i];
            }
            sigma = 1.0/sigma;
//             printf("sigma: %.4f\n", sigma);
            
            for (integer i = 0; i < n * n; i++) eta[i] = 0;
            
            
            for (integer i = 0; i < num; i++)
            {
                // Compute inv(LplusX) -> iLplusX
                dcopy_(&length, I, &GLOBAL::IONE, iLplusX + n * n * i, &GLOBAL::IONE);
                dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, LLplusX + n * n * i, &N, iLplusX + n * n * i, &N, &info);
                
                // compute inv(LplusX) = inv(LplusX)^T * inv(LplusX)
                dcopy_(&length, iLplusX + n * n * i, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
                dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, tmp, &N, &GLOBAL::DZERO, iLplusX + n * n * i, &N);
                
                // add scaled inv(LplusX) to eta
                scale = 0.5/coef[i];
                dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &scale, iLplusX + n * n * i, &N, I, &N, &GLOBAL::DONE, eta, &N);
                
            }
            
            // compute X1 * eta -> tmp
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, X1, &N, eta, &N, &GLOBAL::DZERO, tmp, &N);
            
            // compute tmp * X1 -> eta
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, X1, &N, &GLOBAL::DZERO, eta, &N);
            
            // compute eta - coef * X -> tmp
            dcopy_(&length, eta, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            dcopy_(&length, X1, &GLOBAL::IONE, eta, &GLOBAL::IONE);
            double negnum = 0.0;
            for (integer i = 0; i < num; i++) negnum = negnum - 0.5/coef[i];
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, I, &N, &negnum, eta, &N);            
           
            sigma = 2.0 * sigma/(1.0 - alpha);
            dscal_(&length, &sigma, eta, &GLOBAL::IONE);
           
            
        }
        else
        {
            // Compute inv(LplusX) -> iLplusX
            dcopy_(&length, I, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, LLplusX + n * n * idx, &N, tmp, &N, &info);
            
            // compute inv(LplusX) = inv(LplusX)^T * inv(LplusX)
            dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, tmp, &N, &GLOBAL::DZERO, eta, &N);
            
            // compute X1 * eta -> tmp
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, X1, &N, eta, &N, &GLOBAL::DZERO, tmp, &N);
            
            // compute tmp * X1 -> eta
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, X1, &N, &GLOBAL::DZERO, eta, &N);
            
            // compute eta - X -> eta: -X + eta -> eta
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &negone, X1, &N, I, &N, &GLOBAL::DONE, eta, &N);
            
            scale = 2.0/(1.0 - alpha);
            dscal_(&length, &scale, eta, &GLOBAL::IONE);
        }
        
        
        funsSeries[iter] = 0;
        for (integer i = 0; i < num; i++) funsSeries[iter] = funsSeries[iter] + coef[i];
        
        
//      eta -> -eta
        dscal_(&length, &negone, eta, &GLOBAL::IONE);
        
        
//          --------------  update next iterate  --------------
        dscal_(&length, &stepsize, eta, &GLOBAL::IONE);
        
        
//      Lx^{-1} * eta -> eta
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, L, &N, eta, &N, &info);
//      eta -> eta^T
        dcopy_(&length, eta, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, I, &N, &GLOBAL::DZERO, eta, &N);
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, L, &N, eta, &N, &info);
        
        dcopy_(&length, eta, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, I, &N, &GLOBAL::DZERO, eta, &N);
        
        dcopy_(&length, eta, &GLOBAL::IONE, X2, &GLOBAL::IONE);
        
        Matrix MMt(X2, n, n);
        Matrix::ExpSymmetricM(GLOBAL::L, MMt, MMt);
        
        dcopy_(&length, X2, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, L, &N, tmp, &N, &GLOBAL::DZERO, X2, &N);
        
        dcopy_(&length, X2, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, tmp, &N, L, &N, &GLOBAL::DZERO, X2, &N);
        
        //         update X1 <- X2
        dcopy_(&length, X2, &GLOBAL::IONE, X1, &GLOBAL::IONE);
//         ForDebug::Print("X2", X2, n, n);
        
        
        if(Debug >= 3 && iter > 0)
        {
            disSeries[iter] = distwo(X2, Xtrue, n);
            gradSeries[iter] = 0.0;
        }
        
        timeSeries[iter] = static_cast<double>(getTickCount() - starttime) / CLK_PS; // record the time for each iteration
        
        // check stopping criteria
        if(iter == 0)
        {
            err_g = std::abs(funsSeries[0]);
        }
        else
        {
            err_g = std::abs(funsSeries[iter] - funsSeries[iter - 1])/funsSeries[0];
        }
        enditer = iter;
        
        
//         if(err_g < TOL)
//         {
//             enditer = iter;
//             printf("error: %.4f\n", err_g);
//             printf("Tolerance reached!");
//             break;
//         }
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
    plhs[5] = mxCreateDoubleMatrix(maxiter, 1, mxREAL); // function
    plhs[6] = mxCreateDoubleMatrix(maxiter, 1, mxREAL); // grad
    plhs[7] = mxCreateDoubleScalar(err_g);
    
    double *plhstime = mxGetPr(plhs[0]);
    double *plhsdis = mxGetPr(plhs[4]);
    double *plhsfun = mxGetPr(plhs[5]);
    double *plhsgrad = mxGetPr(plhs[6]);
    
    for(integer i = 0; i < maxiter; i++)
    {
        plhstime[i] = timeSeries[i];
        plhsdis[i] = disSeries[i];
        plhsfun[i] = funsSeries[i];
        plhsgrad[i] = gradSeries[i];
//         printf("there: %.4e, %.4e, %.4e, %.4e\n", timeSeries[i], disSeries[i], stepsizeSeries[i], truestepsizeSeries[i]);
    }
    
    delete[] timeSeries;
    delete[] disSeries;
    delete[] gradSeries;
    delete[] funsSeries;
    
    delete[] coef;
    delete[] eta;
    delete[] L;
    delete[] tmp;
    delete[] X1;
    delete[] X2;
    delete[] I;
    delete[] CLplusX;
    delete[] LLplusX;
    delete[] iLplusX;
    
    
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







