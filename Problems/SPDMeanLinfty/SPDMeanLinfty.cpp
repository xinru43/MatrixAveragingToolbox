
#include "SPDMeanLinfty.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDMeanLinfty::SPDMeanLinfty(double *inLs, integer inn, integer innum)
	{
		Ls = inLs;
		n = inn;
		num = innum;
	};

	SPDMeanLinfty::~SPDMeanLinfty(void)
	{
	};

	double SPDMeanLinfty::f(Variable *x) const
	{
		SPDManifold *Mani = dynamic_cast<SPDManifold *> (Domain);
		if (!x->TempDataExist("L"))
		{
			Mani->CholeskyRepresentation(x);
		}
		const SharedSpace *SharedL = x->ObtainReadTempData("L");
		Variable *LElement = SharedL->GetSharedElement();
		const double *L = LElement->ObtainReadData();

		SharedSpace *SharedlogLXL = new SharedSpace(3, n, n, num);
		double *logLXL = SharedlogLXL->ObtainWriteEntireData();
        
		double *Ltmp = new double[n * n];
		integer length = n * n, N = n, info;
        
        
        double *coef = new double[num];
        
        SharedSpace *Sharedmaxval = new SharedSpace(1, 2); 
        double *maxval = Sharedmaxval->ObtainWriteEntireData();

        
        
		for (integer i = 0; i < num; i++)
		{
			dcopy_(&length, const_cast<double *> (L), &GLOBAL::IONE, Ltmp, &GLOBAL::IONE);
			/*Solve the linear system Li X = Lx, i.e., X = Li^{-1} Lx. The solution X is stored in LiiLx.
			Note that Li is a lower triangular matrix.
			Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
			dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, Ls + n * n * i, &N, Ltmp, &N, &info);
			if (info != 0)
			{
				printf("The cholesky decompsotion in SPDMean::f failed with info:%d!\n", info);
			}
			dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, Ltmp, &N, Ltmp, &N, &GLOBAL::DZERO, logLXL + n * n * i, &N);
			Matrix MMt(logLXL + n * n * i, n, n);
			Matrix::LogSymmetricM(GLOBAL::L, MMt, MMt);
            
            coef[i] = dnrm2_(&length, logLXL + n * n * i, &GLOBAL::IONE);
		}
        
        
        double result = coef[0];
        maxval[1] = 0;
        for (integer i = 1; i < num; i++)
        {
            if(coef[i] > result)
            {
                result = coef[i];
                maxval[1] = static_cast<double>(i);
//                 printf("I am the true index: %d, I am not: %d\n", i, static_cast<integer>(maxval[1]));
            }
        }
        
        
        
        maxval[0] = result;
//         printf("I am the true result: %.2e, I am not: %.2e\n", result, maxval[0]);

        delete[] Ltmp;
        delete[] coef;
		x->AddToTempData("logLXL", SharedlogLXL);
        x->AddToTempData("maxval", Sharedmaxval);

		return result;
  
// -------------------------------------------------------------------------
//         SPDManifold *Mani = dynamic_cast<SPDManifold *> (Domain);
// 		if (!x->TempDataExist("L"))
// 		{
// 			Mani->CholeskyRepresentation(x);
// 		}
// 		const SharedSpace *SharedL = x->ObtainReadTempData("L");
// 		Variable *LElement = SharedL->GetSharedElement();
// 		const double *L = LElement->ObtainReadData();
// 
// 		SharedSpace *SharedlogLXL = new SharedSpace(3, n, n, num);
// 		double *logLXL = SharedlogLXL->ObtainWriteEntireData();
// 		double *Ltmp = new double[n * n];
// 		integer length = n * n, N = n, info;
//         
//         
//         double *eigenvalues = new double[n * num];
//         double *eigenvectors = new double[2 * n * n];
//         double *eigenvectorsD = eigenvectors + n * n;
//         
//         double a, result;
//         SharedSpace *Sharedcoef = new SharedSpace(3, num, 1, 1); 
//         double *coef = Sharedcoef->ObtainWriteEntireData();
//         
//         
// 		for (integer i = 0; i < num; i++)
// 		{
// 			dcopy_(&length, const_cast<double *> (L), &GLOBAL::IONE, Ltmp, &GLOBAL::IONE);
// 			/*Solve the linear system Li X = Lx, i.e., X = Li^{-1} Lx. The solution X is stored in LiiLx.
// 			Note that Li is a lower triangular matrix.
// 			Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
// 			dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, Ls + n * n * i, &N, Ltmp, &N, &info);
// 			if (info != 0)
// 			{
// 				std::cout << "The cholesky decompsotion in SPDMean::f failed with info:" << info << "!" << std::endl;
// 			}
// 			dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, Ltmp, &N, Ltmp, &N, &GLOBAL::DZERO, logLXL + n * n * i, &N);
//             
//              Matrix E(eigenvalues + n * i, n, 1), V(eigenvectors, n, n), VD(eigenvectorsD, n, n);
//              Matrix MMt(logLXL + n * n * i, n, n);
//              Matrix::EigenSymmetricM(GLOBAL::L, MMt, E, V);
// 
//              dcopy_(&length, eigenvectors, &GLOBAL::IONE, VD.matrix, &GLOBAL::IONE);
//              for (integer j = 0; j < n; j++)
//              {
//                  a = log(eigenvalues[j + n * i]);
//                  dscal_(&N, &a, eigenvectors + j * n, &GLOBAL::IONE);
//              }
//              
//              Matrix MlogLXL(logLXL + n * n * i, n, n);
//              Matrix::DGEMM(GLOBAL::DONE, V, false, VD, true, GLOBAL::DZERO, MlogLXL);
//              
//              coef[i] = 0.0;
//              for(integer j = 0; j < N; j++) {coef[i] = coef[i] + log(eigenvalues[j + i * N]) * log(eigenvalues[j + i * N]);}
//              coef[i] = std::sqrt(coef[i]);
// 		}
//         
//         result = coef[0];
//         for (integer i = 1; i < num; i++)
//         {
//             if (coef[i] > result)
//             {
//                 result = coef[i];
//             }
//         }
//         
// 
// 		delete[] Ltmp;
//         delete[] eigenvalues;
//         delete[] eigenvectors;
// 		x->AddToTempData("logLXL", SharedlogLXL);
//         x->AddToTempData("coef", Sharedcoef);        
// 
// 		return result;
	};
    
    

	void SPDMeanLinfty::RieGrad(Variable *x, Vector *gf) const
	{
		const SharedSpace *SharedlogLXL = x->ObtainReadTempData("logLXL");
		const double *logLXL = SharedlogLXL->ObtainReadData();
		double *gfVT = gf->ObtainWriteEntireData();
		for (integer i = 0; i < n * n; i++)
			gfVT[i] = 0;
		const double *xM = x->ObtainReadData();

		integer N = n, info;

        const SharedSpace *Sharedmaxval = x->ObtainReadTempData("maxval");
        const double *maxval  = Sharedmaxval->ObtainReadData();
        
        double *tmp = new double[n * n];

        integer idx = static_cast<integer>(maxval[1]);
        
        dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, const_cast<double *> (logLXL + n * n * idx), &N,
                Ls + n * n * idx, &N, &GLOBAL::DZERO, tmp, &N);
        
        /*Solve the linear system Li^T X = tmp, i.e., X = Li^{-T} log(Li^{-1} X Li^{-T}) Li^T. The solution X is stored in tmp.
         * Note that Li is a lower triangular matrix.
         * Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
        dtrtrs_(GLOBAL::L, GLOBAL::T, GLOBAL::N, &N, &N, Ls + n * n * idx, &N, tmp, &N, &info);
        if (info != 0)
        {
            printf("The cholesky decompsotion in SPDMean::RieGrad failed with info:%d!\n", info);
        }
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *> (xM), &N, tmp, &N, &GLOBAL::DONE, gfVT, &N);

        
        double coefa = 1.0/maxval[0];
        integer length = n * n;
        dscal_(&length, &coefa, gfVT, &GLOBAL::IONE);
        delete[] tmp;
    };
    
    

	void SPDMeanLinfty::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const
	{
		std::cout << "warning: SPDMean::RieHessianEta has not been implemented!" << std::endl;
		etax->CopyTo(xix);
	};
}; /*end of ROPTLIB namespace*/
