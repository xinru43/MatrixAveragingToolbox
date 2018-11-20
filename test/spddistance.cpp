
/*Help to debug the code*/
// #include "ForDebug.h"
// #include <iostream>
// #include <ctime>
// #include "Timer.h"
// #include "def.h"
// #include "MyMatrix.h"
// #include "SPDMean.h"
// #include "SolversLS.h"

#include "spddistance.h"

// using namespace ROPTLIB;

// void geodesic_t(double *A, double *B, integer n, double t, double *output);

#if !defined(MATLAB_MEX_FILE) && defined(TESTSPDDISTANCE)
int main()
{
    testspddistance();
#ifdef _WIN64
#ifdef _DEBUG
    _CrtDumpMemoryLeaks();
#endif
#endif
    return 0;
}
#endif





void testspddistance(void)
{
    printf("============= RSD: Using C++ =============\n");
// 	// choose a random seed
// 	unsigned tt = (unsigned)time(NULL);
// 	tt = 0;
// //	genrandseed(tt);
//
// 	/*Randomly generate a point on the SPD manifold*/
// 	integer n = 3, num = 3;
//     double *initial_X = new double[n * n];
// 	for (integer i = 0; i < n; i++)
// 	{
// 		for (integer j = 0; j < n; j++)
// 		{
// 			initial_X[i + j * n] = 0;
// 		}
// 		initial_X[i + i * n] = 1;
// 	}
};




#ifdef MATLAB_MEX_FILE

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
//     printf("============= RSD: Using mex =============\n");
    // deal with input
    double *X, *Y;
    integer n;
    double alpha;
    X = mxGetPr(prhs[0]);
    Y = mxGetPr(prhs[1]);
    n = static_cast<integer>(mxGetScalar(prhs[2]));
    char dis_type[20] = "";
    mxGetString(prhs[3], dis_type, 20);
    alpha = static_cast<double>(mxGetScalar(prhs[4]));
    

    unsigned long starttime = getTickCount();
    double result, ComTime;
    integer N = n, length = n * n;
    std::string stddis_type = dis_type;
    
//     make identity matrix
    double *I = new double[length];
    for(integer i = 0; i < length; i++) I[i] = 0.0;
    for(integer i = 0; i < n; i++) I[i*n+i] = 1.0;
    
    
//------------------------------- Euclidean distance -------------------------------
    if (stddis_type == "euclid")
    {
        double diff;
        result = 0.0;
        for (integer i = 0; i < length; i++)
        {
//             error = const_cast<double *>(x1M + i) - const_cast<double *>(x2M + i);
            diff = X[i] - Y[i];
            result = result + diff * diff;
        }
        
        result = sqrt(result);
    }
    
    
//------------------------------- Log-Euclidean distance -------------------------------
    if (stddis_type == "logeuclid")
    {
        double *logX = new double[n * n];
        double *logY = new double[n * n];
        
//         compute log(X)
        dcopy_(&length, X, &GLOBAL::IONE, logX, &GLOBAL::IONE);
        Matrix MMx(logX, n, n);
        Matrix::LogSymmetricM(GLOBAL::L, MMx, MMx);
//         compute log(Y)
        dcopy_(&length, Y, &GLOBAL::IONE, logY, &GLOBAL::IONE);
        Matrix MMy(logY, n, n);
        Matrix::LogSymmetricM(GLOBAL::L, MMy, MMy);
        
        result = 0.0;
        double diff;
        for (integer i = 0; i < length; i++)
        {
            diff = logX[i] - logY[i];
            result = result + diff * diff;
        }
        result = sqrt(result);
        
        delete[] logX;
        delete[] logY;
    }
    
    
    
    
    //------------------------------- Riemannian distance -------------------------------
    if (stddis_type == "riemann")
    {
        // dis(x, y) = |log(Lx^{-1} y Lx^{-T})|_F
        
        // Cholesky factorization of X
        integer info;
        double *Lx = new double[n * n];
        double *Ly = new double[n * n];
        double *LxLy = new double[n * n];
        dcopy_(&length, X, &GLOBAL::IONE, Lx, &GLOBAL::IONE);
        dcopy_(&length, Y, &GLOBAL::IONE, Ly, &GLOBAL::IONE);
        
        dpotrf_(GLOBAL::L, &N, Lx, &N, &info);
        dpotrf_(GLOBAL::L, &N, Ly, &N, &info);
        for (integer j = 0; j < n; j++)
        {
            for (integer k = j + 1; k < n; k++)
            {
                Lx[j + k * n] = 0;
                Ly[j + k * n] = 0;
            }
        }
        
        dcopy_(&length, Ly, &GLOBAL::IONE, LxLy, &GLOBAL::IONE);
        /*Solve the linear system Lx X = Ly, i.e., X = Lx^{-1}Ly. The solution X is stored in LxLy*/
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, Lx, &N, LxLy, &N, &info);
        
        /* compute temp = LxLy*(LxLy)^T  */
        double *temp = new double[n * n];
        dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, LxLy, &N, LxLy, &N, &GLOBAL::DZERO, temp, &N);
        
        double *eigenvalues = new double[n + n * n];
        double *eigenvectors = eigenvalues + n;
        Matrix E(eigenvalues, n, 1), V(eigenvectors, n, n);
        Matrix MMt(temp, n, n);
        Matrix::EigenSymmetricM(GLOBAL::L, MMt, E, V);
        
        result = 0.0;
        for (integer i = 0; i < N; i++)
        {
            result = result + log(eigenvalues[i])*log(eigenvalues[i]);
        }
        
        result = sqrt(result);
        
        delete[] eigenvalues;
        delete[] temp;
        delete[] LxLy;
        delete[] Lx;
        delete[] Ly;
    }
    
    
    
    
    
    
    //------------------------------- Alpha divergence -------------------------------
    if (stddis_type == "alpha")
    {
        double a1 = (1.0 - alpha)/2.0, a2 = (1.0 + alpha)/2.0;
        double *CXplusY = new double[n * n];
        integer info;
        dcopy_(&length, Y, &GLOBAL::IONE, CXplusY, &GLOBAL::IONE);
        // compute a1 * X + a2 * Y -> CXplusY
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &a1, X, &N, I, &N, &a2, CXplusY, &N);
        
        // Cholesky factorization of CXplusY
        dpotrf_(GLOBAL::L, &N, CXplusY, &N, &info);
        
        // compute the determinant
        result = 1.0;
        for (integer j = 0; j < n; j++)
        {
            result = result * CXplusY[j * n + j];
        }
        result = log(result * result);
        
        // compute the determinant of X: product of diagonal elements
        double detX = 1.0, detY = 1.0;
        double *Lx = new double[n * n];
        double *Ly = new double[n * n];
        dcopy_(&length, X, &GLOBAL::IONE, Lx, &GLOBAL::IONE);
        dpotrf_(GLOBAL::L, &N, Lx, &N, &info);
        dcopy_(&length, Y, &GLOBAL::IONE, Ly, &GLOBAL::IONE);
        dpotrf_(GLOBAL::L, &N, Ly, &N, &info);
        
        for (integer j = 0; j < n; j++)
        {
            detX = detX * Lx[j * n + j];
            detY = detY * Ly[j * n + j];
        }
//         printf("S1: detX: %.2f\n", detX);
//         printf("S1: detY: %.2f\n", detY);
        detX = a1 * log(std::abs(detX * detX));
        detY = a2 * log(std::abs(detY * detY));
//         printf("S2: detX: %.2f\n", detX);
//         printf("S2: detY: %.2f\n", detY);
        
        result = result - detX - detY;
        result = (result * 4.0)/(1 - alpha * alpha);
        result = sqrt(std::abs(result));
        
        delete[] CXplusY;
        delete[] Lx;
        delete[] Ly;
    }
    
    
    
        //------------------------------- Symmetric Alpha divergence -------------------------------
    if (stddis_type == "symmalpha")
    {
        double a1 = (1.0 - alpha)/2.0, a2 = (1.0 + alpha)/2.0;
        double *CXplusY = new double[n * n];
        double *CYplusX = new double[n * n];
        integer info;
        dcopy_(&length, Y, &GLOBAL::IONE, CXplusY, &GLOBAL::IONE);
        // compute a1 * X + a2 * Y -> CXplusY
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &a1, X, &N, I, &N, &a2, CXplusY, &N);
        dpotrf_(GLOBAL::L, &N, CXplusY, &N, &info);  // Cholesky factorization of CXplusY
        
        // compute the determinant
        double result1 = 1.0;
        for (integer j = 0; j < n; j++)
        {
            result1 = result1 * CXplusY[j * n + j];
        }
        result1 = log(result1 * result1);
    
        dcopy_(&length, Y, &GLOBAL::IONE, CYplusX, &GLOBAL::IONE);
        // compute a2 * X + a1 * Y -> CYplusX
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &a2, X, &N, I, &N, &a1, CYplusX, &N);
        dpotrf_(GLOBAL::L, &N, CYplusX, &N, &info);
        
        // compute the determinant
        double result2 = 1.0;
        for (integer j = 0; j < n; j++)
        {
            result2 = result2 * CYplusX[j * n + j];
        }
        result2 = log(result2 * result2);
    
        // compute the determinant of X & Y: product of diagonal elements
        double detX = 1.0, detY = 1.0;
        double *Lx = new double[2 * n * n];
        double *Ly = Lx + n * n;
        dcopy_(&length, X, &GLOBAL::IONE, Lx, &GLOBAL::IONE);
        dpotrf_(GLOBAL::L, &N, Lx, &N, &info);
        dcopy_(&length, Y, &GLOBAL::IONE, Ly, &GLOBAL::IONE);
        dpotrf_(GLOBAL::L, &N, Ly, &N, &info);
        
        for (integer j = 0; j < n; j++)
        {
            detX = detX * Lx[j * n + j];
            detY = detY * Ly[j * n + j];
        }
        detX = log(detX * detX);
        detY = log(detY * detY);
        
        result = result1 + result2 - detX - detY;
        result = sqrt(std::abs(result));
        
        delete[] CXplusY;
        delete[] CYplusX;
        delete[] Lx;
    }
    
    
    
        //------------------------------- Symmetric LD Bregman divergence/Jeffrey -------------------------------
    if (stddis_type == "symmldbreg")
    {
        integer info;
        double *Lx = new double[n * n];
        double *Ly = new double[n * n];
        double *tmp = new double[n * n];
        dcopy_(&length, X, &GLOBAL::IONE, Lx, &GLOBAL::IONE);
        dcopy_(&length, Y, &GLOBAL::IONE, Ly, &GLOBAL::IONE);
        
        dpotrf_(GLOBAL::L, &N, Lx, &N, &info);
        dpotrf_(GLOBAL::L, &N, Ly, &N, &info);
        for (integer j = 0; j < n; j++)
        {
            for (integer k = j + 1; k < n; k++)
            {
                Lx[j + k * n] = 0;
                Ly[j + k * n] = 0;
            }
        }
    
        dcopy_(&length, X, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        /*Solve the linear system Ly * tmp = X, i.e., tmp = Ly^{-1}X. The solution is stored in tmp*/
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, Ly, &N, tmp, &N, &info);
        dtrtrs_(GLOBAL::L, GLOBAL::T, GLOBAL::N, &N, &N, Ly, &N, tmp, &N, &info); 
//         ForDebug::Print("tmp", tmp, n, n);
        
        result = 0.0;
        for (integer i = 0; i < n; i++)
        {
            result = result + tmp[i * n + i];
        }
//         printf("1: %.2f\n", result);

        dcopy_(&length, Y, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        /*Solve the linear system Ly * tmp = X, i.e., tmp = Ly^{-1}X. The solution is stored in tmp*/
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, Lx, &N, tmp, &N, &info);
        dtrtrs_(GLOBAL::L, GLOBAL::T, GLOBAL::N, &N, &N, Lx, &N, tmp, &N, &info); 
//         ForDebug::Print("tmp", tmp, n, n);

        for (integer i = 0; i < n; i++)
        {
            result = result + tmp[i * n + i];
        }
//         printf("2: %.2f\n", result);
        
        result = result - 2.0 * n;
//         printf("3: %.2f\n", result);
        result = sqrt(std::abs(result));
//         printf("4: %.2f\n", result);
        delete[] tmp;
        delete[] Lx;
        delete[] Ly;
    }
    
    
    
    
    
    ComTime = static_cast<double>(getTickCount() - starttime) / CLK_PS;
    
    plhs[0] = mxCreateDoubleScalar(result);
    plhs[1] = mxCreateDoubleScalar(ComTime);
    
    delete[] I;
    
    if(isnan(result)) 
    {
        result = 0.0;        
    }
    
    return;
    
}

#endif













