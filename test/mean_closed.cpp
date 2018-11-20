
/*Help to debug the code*/
// #include "ForDebug.h"
// #include <iostream>
// #include <ctime>
// #include "Timer.h"
// #include "def.h"
// #include "MyMatrix.h"
// #include "SPDMean.h"
// #include "SolversLS.h"

#include "mean_closed.h"

// using namespace ROPTLIB;

void geodesic_t(double *A, double *B, integer n, double t, double *output);

#if !defined(MATLAB_MEX_FILE) && defined(TESTSPDDLONEPARAMFP)
int main()
{
     testmeanclosed();
#ifdef _WIN64
#ifdef _DEBUG
    _CrtDumpMemoryLeaks();
#endif
#endif
    return 0;
}
#endif





void testmeanclosed(void)
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
    double *As;
    integer n, num;
    double alpha;
    As = mxGetPr(prhs[0]);
    n = static_cast<integer>(mxGetScalar(prhs[1]));
    num = static_cast<integer>(mxGetScalar(prhs[2]));
    
    char mean_type[20] = "";
	mxGetString(prhs[3], mean_type, 20);
    
   
    unsigned long starttime = getTickCount();
    double ComTime;
    integer length = n * n;
    std::string stdmean_type = mean_type;
    double *output = new double[length];
    
//     make identity matrix
    double *I = new double[length];
    for(integer i = 0; i < length; i++) I[i] = 0.0;
    for(integer i = 0; i < n; i++) I[i*n+i] = 1.0;
    

//------------------------------- Euclidean mean -------------------------------  
    if (stdmean_type == "euclid")
    {
        for(integer j = 0; j < length; j++) output[j] = As[j];
        
        for(integer i = 1; i < num; i++)
        {
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, As + n * n * i, &n, I, &n, &GLOBAL::DONE, output, &n);
        }
        double a = 1.0/num;
        dscal_(&length, &a, output, &GLOBAL::IONE); 
    }
    
    
//------------------------------- Log-Euclidean mean ------------------------------- 
    if (stdmean_type == "logeuclid")
    {
        double *logAs = new double[n * n];
        for(integer j = 0; j < length; j++) output[j] = 0.0;
        
        for(integer i = 0; i < num; i++)
        {
            dcopy_(&length, As + n * n * i, &GLOBAL::IONE, logAs, &GLOBAL::IONE);
            Matrix MMt(logAs, n, n);
            Matrix::LogSymmetricM(GLOBAL::L, MMt, MMt);
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, logAs, &n, I, &n, &GLOBAL::DONE, output, &n);
        }
        double a = 1.0/num;
        dscal_(&length, &a, output, &GLOBAL::IONE); 
        Matrix MMt(output, n, n);
        Matrix::ExpSymmetricM(GLOBAL::L, MMt, MMt);
        delete[] logAs;
    }
    
    
//------------------------------- Inductive mean -------------------------------
    if (stdmean_type == "inductive")
    {
        double *tmp = new double[n * n];
        double t = 1.0/2;
        geodesic_t(As, As + n * n, n, t, output);
        for(integer i = 2; i < num; i++)
        {
            t = 1.0/(i + 1);
            geodesic_t(output, As + n * n * i, n, t, tmp);
            dcopy_(&length, tmp, &GLOBAL::IONE, output, &GLOBAL::IONE);
        }
        delete[] tmp;
        
    }
    
//------------------------------- Arithmetic-Harmonic mean -------------------------------
    if (stdmean_type == "symmldbreg")
    {
        
        // compute the arithmetic mean
        double *Ar = new double[n * n];
        double *Ha = new double[n * n];
        
        dcopy_(&length, As, &GLOBAL::IONE, Ar, &GLOBAL::IONE);
        for(integer i = 1; i < num; i++)
        {
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, As + n * n * i, &n, I, &n, &GLOBAL::DONE, Ar, &n);
        }
        double a = 1.0/num;
        dscal_(&length, &a, Ar, &GLOBAL::IONE);
        
        
        // compute the harmonic mean
        double *Ls = new double[n * n * num];
        double *tmp = new double[n * n];
        integer info;
        integer Nnum = n * n * num;
        for(integer i = 0; i < length; i++) Ha[i] = 0.0;
        
        dcopy_(&Nnum, As, &GLOBAL::IONE, Ls, &GLOBAL::IONE);
        for(integer i = 0; i < num; i++)
        {
            // Cholesky factorization of CLplusX
            dpotrf_(GLOBAL::L, &n, Ls + n * n * i, &n, &info);
            for (integer j = 0; j < n; j++)
                for (integer k = j + 1; k < n; k++)
                    Ls[j + k * n + i * n * n] = 0.0;
            
            // Compute inv(Ls) -> tmp
            dcopy_(&length, I, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &n, &n, Ls + n * n * i, &n, tmp, &n, &info);
            
            // compute inv(Ls) = inv(Ls)^T * inv(Ls): tmp^T * tmp -> Ls
            dgemm_(GLOBAL::T, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, tmp, &n, tmp, &n, &GLOBAL::DZERO, Ls + n * n * i, &n);
            
            // add inv(Ls) to Ha
            dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, Ls + n * n * i, &n, I, &n, &GLOBAL::DONE, Ha, &n);
        }
        
        // Cholesky factorization of Ha
        dpotrf_(GLOBAL::L, &n, Ha, &n, &info);
        for (integer j = 0; j < n; j++)
            for (integer k = j + 1; k < n; k++)
                Ha[j + k * n] = 0.0;
        
        // compute inv(L_Ha)
        dcopy_(&length, I, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &n, &n, Ha, &n, tmp, &n, &info);
        
        // compute inv(Ha) = inv(L_Ha)^T * inv(L_Ha)
        dgemm_(GLOBAL::T, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, tmp, &n, tmp, &n, &GLOBAL::DZERO, Ha, &n);
        double dnum = 1.0 * num;
        dscal_(&length, &dnum, Ha, &GLOBAL::IONE);        
        
        // compute the geometric mean of Ar and Ha
        geodesic_t(Ha, Ar, n, 0.5, output); 
        
        delete[] Ar;
        delete[] Ha;
        delete[] tmp;
        delete[] Ls;
        
    }    
    

    ComTime = static_cast<double>(getTickCount() - starttime) / CLK_PS;
    
    
    // create output
    plhs[0] = mxCreateDoubleMatrix(length, 1, mxREAL); // X
    double *covmean = mxGetPr(plhs[0]);
    dcopy_(&length, output, &GLOBAL::IONE, covmean, &GLOBAL::IONE);
    plhs[1] = mxCreateDoubleScalar(ComTime); //  total computational time
    
    delete[] output;
    delete[] I;

    return;
    
}

#endif




void geodesic_t(double *A, double *B, integer n, double t, double *output)
{
   double *a = new double[2 * n * n];
   double *b = a + n * n;
   integer info;
   
    // Cholesky factorization of A, B
   integer length = n * n;
   dcopy_(&length, A, &GLOBAL::IONE, a, &GLOBAL::IONE);
   dcopy_(&length, B, &GLOBAL::IONE, b, &GLOBAL::IONE);

   
   dpotrf_(GLOBAL::L, &n, a, &n, &info);
   for (integer j = 0; j < n; j++)
       for (integer k = j + 1; k < n; k++)
           a[j + k * n] = 0;
   
   dpotrf_(GLOBAL::L, &n, b, &n, &info);
   for (integer j = 0; j < n; j++)
       for (integer k = j + 1; k < n; k++)
           b[j + k * n] = 0;
    
//    solve a * tmp = b: tmp <- a^(-1) * b
   double *tmp = new double[n * n];
   dcopy_(&length, b, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
   dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &n, &n, a, &n, tmp, &n, &info);
   

//    compute tmp * tmp' -> C
   double *C = new double[n * n];
   dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, tmp, &n, tmp, &n, &GLOBAL::DZERO, C, &n);
   
//    eigenvalue decomposition of C
		double *eigenvalues = new double[n + 2 * n * n];
		double *eigenvectors = eigenvalues + n;
		double *eigenvectorsD = eigenvectors + n * n;
        Matrix S(C, n, n), Ct(output, n, n);
		Matrix E(eigenvalues, n, 1), V(eigenvectors, n, n), VD(eigenvectorsD, n, n);
		Matrix::EigenSymmetricM(GLOBAL::L, S, E, V);
		dcopy_(&length, eigenvectors, &GLOBAL::IONE, VD.matrix, &GLOBAL::IONE);
        
        double coef;
		for (integer i = 0; i < n; i++)
		{
			coef = std::pow(eigenvalues[i], t);
			dscal_(&n, &coef, eigenvectors + i * n, &GLOBAL::IONE);
		}
		Matrix::DGEMM(GLOBAL::DONE, V, false, VD, true, GLOBAL::DZERO, Ct);
        

        dgemm_(GLOBAL::N, GLOBAL::N, &n, &n, &n, &GLOBAL::DONE, a, &n, output, &n, &GLOBAL::DZERO, tmp, &n);
        dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, tmp, &n, a, &n, &GLOBAL::DZERO, output, &n);
    
         
        delete[] a;
        delete[] C;
		delete[] eigenvalues;
        delete[] tmp;
};












