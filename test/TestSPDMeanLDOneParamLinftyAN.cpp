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

#include "TestSPDMeanLDOneParamLinftyAN.h"

using namespace ROPTLIB;

double dis(double *X, integer n);
double distwo(double *X, double *Y, integer n);

#if !defined(MATLAB_MEX_FILE) && defined(TESTSPDMEANLDONEPARAMLINFTYAN)

int main()
{
    testSPDMeanLDOneParamLinftyAN();
    
#ifdef _WIN64
#ifdef _DEBUG
    _CrtDumpMemoryLeaks();
#endif
#endif
    return 0;
}

#endif


void testSPDMeanLDOneParamLinftyAN(void)
{
    printf("============= LogDet alpha-divergence minimax center computation by AN method =============\n");
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
    
    double maxval = -1.0, negone = -1.0, t;
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
        
        
        maxval = coef[0];
        idx = 0;
        for (integer i = 1; i < num; i++)
        {
            if(coef[i] > maxval)
            {
                maxval = coef[i];
                idx = i;
            }
        }
        
        funsSeries[iter] = maxval;
        
//          --------------  compute eta  --------------
//          logLXL * L^T -> tmp
        dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, logLXL + n * n * idx, &N, Ls + n * n * idx, &N, &GLOBAL::DZERO, tmp, &N);
        
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
    
    printf("============= LogDet alpha-divergence minimax center computation by AN method =============\n");
    // deal with input
    double *Ls, *LDAs, *initial_X, *Xtrue;
    integer n, num, maxiter, Debug;
    double TOL, alpha;
    Ls = mxGetPr(prhs[0]);
    LDAs = mxGetPr(prhs[1]);
    alpha = static_cast<double>(mxGetScalar(prhs[2]));
    initial_X = mxGetPr(prhs[3]);
    Xtrue = mxGetPr(prhs[4]);
    n = static_cast<integer>(mxGetScalar(prhs[5]));
    num = static_cast<integer>(mxGetScalar(prhs[6]));
    maxiter = static_cast<integer>(mxGetScalar(prhs[7]));
    TOL = static_cast<double>(mxGetScalar(prhs[8]));
    Debug = static_cast<integer>(mxGetScalar(prhs[9]));
    
    
//     time start
    unsigned long starttime = getTickCount();
    CheckMemoryDeleted = new std::map<integer *, integer>;
    
    integer length = n * n, N = n, Nnum = n * n * num, info, enditer, idx;
    double maxval = -1.0, negone = -1.0, t;
    double ComTime, err_g;
    double a1 = (1.0 - alpha)/2.0, a2 = (1.0 + alpha)/2.0;
    
    double *X1 = new double[n * n];
    dcopy_(&length, initial_X, &GLOBAL::IONE, X1, &GLOBAL::IONE);
    
    double *timeSeries = new double[maxiter];
    double *funsSeries = new double[maxiter];
    double *gradSeries = new double[maxiter];
    double *disSeries = new double[maxiter];
    
    double *logLXL = new double[n * n];
    double *CLplusX = new double[n * n * num];
    double *coef = new double[num];
    double *eta = new double[n * n];
    double *Lx = new double[n * n];
    double *tmp = new double[n * n];
    double *X2 = new double[n * n];
    
    double *In = new double[n * n]; // make identity matrix
    for (integer i = 0; i < length; i++) In[i] = 0.0;
    for (integer i = 0; i < n; i++) In[i*n+i] = 1.0;
    
    
    if(Debug >= 3)
    {
        // compute distance between current iterate and identity matrix
        disSeries[0] = distwo(X1, Xtrue, n);
        gradSeries[0] = 0.0;
    }
    
    // restore As from cholesky of As
    double *As = new double[n * n * num];
    for (integer i = 0; i < num; i++)
    {
        dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, Ls + n * n * i, &N, Ls + n * n * i, &N, &GLOBAL::DZERO, As + n * n * i, &N);
    }
    
    
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
        dcopy_(&Nnum, As, &GLOBAL::IONE, CLplusX, &GLOBAL::IONE);
        
//   --------------  divergence from A_i to X_k  --------------
        for (integer i = 0; i < num; i++)
        {
            // compute a1 * Ai + a2 * x
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &a2, X1, &N, In, &N, &a1, CLplusX + n * n * i, &N);
            
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
        double detX = 1.0;
        for (integer j = 0; j < n; j++)
        {
            detX = detX * Lx[j * n + j];
        }
        detX = log(detX * detX);
        
        
        for (integer i = 0; i < num; i++)
        {
            coef[i] = coef[i] - a2 * detX - a1 * LDAs[i];
        }
        
        maxval = coef[0];
        idx = 0;
        for (integer i = 1; i < num; i++)
        {
            if(coef[i] > maxval)
            {
                maxval = coef[i];
                idx = i;
            }
        }
        
        funsSeries[iter] = maxval;
        
//          --------------  compute eta  --------------
        dcopy_(&length, Lx, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        /*Solve the linear system Li X = Lx, i.e., X = Li^{-1} Lx. The solution X is stored in LiiLx.
         * Note that Li is a lower triangular matrix.
         * Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, Ls + n * n * idx, &N, tmp, &N, &info);
        if (info != 0)
        {
            printf("The cholesky decompsotion in SPDMean::f failed with info:%d!\n", info);
        }
        dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, tmp, &N, tmp, &N, &GLOBAL::DZERO, logLXL, &N);
        Matrix MMt(logLXL, n, n);
        Matrix::LogSymmetricM(GLOBAL::L, MMt, MMt);
        
//          logLXL * Li^T -> tmp
        dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, logLXL, &N, Ls + n * n * idx, &N, &GLOBAL::DZERO, tmp, &N);
        
//          Li^{-T} * tmp -> tmp, i.e.,  Li^{-T} * logLXL * Li^T -> tmp
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
        
        
//      Lx^{-1} * eta -> eta
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, Lx, &N, eta, &N, &info);
//      eta -> eta^T
        dcopy_(&length, eta, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, In, &N, &GLOBAL::DZERO, eta, &N);
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, Lx, &N, eta, &N, &info);
        
        dcopy_(&length, eta, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, In, &N, &GLOBAL::DZERO, eta, &N);        
        
        dcopy_(&length, eta, &GLOBAL::IONE, X2, &GLOBAL::IONE);
        
        Matrix MMt2(X2, n, n);
        Matrix::ExpSymmetricM(GLOBAL::L, MMt2, MMt2);
        
        dcopy_(&length, X2, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, Lx, &N, tmp, &N, &GLOBAL::DZERO, X2, &N);
        
        dcopy_(&length, X2, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, tmp, &N, Lx, &N, &GLOBAL::DZERO, X2, &N);
        
        //         update X1 <- X2
        dcopy_(&length, X2, &GLOBAL::IONE, X1, &GLOBAL::IONE);        
        
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
    
    delete[] logLXL;
    delete[] CLplusX;
    delete[] As;
    delete[] coef;
    delete[] eta;
    delete[] Lx;
    delete[] tmp;
    delete[] X1;
    delete[] X2;
    delete[] In;
    
    
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







