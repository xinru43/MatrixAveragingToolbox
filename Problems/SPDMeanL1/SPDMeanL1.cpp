
#include "SPDMeanL1.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDMeanL1::SPDMeanL1(double *inLs, integer inn, integer innum)
	{
		Ls = inLs;
		n = inn;
		num = innum;
	};

	SPDMeanL1::~SPDMeanL1(void)
	{
	};

	double SPDMeanL1::f(Variable *x) const
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
		
		integer length = n * n, N = n, info;
        double a, result = 0.0;
        
        double *Ltmp = new double[n * n];
        double *eigenvalues = new double[n * num];
        double *eigenvectors = new double[2 * n * n];
        double *eigenvectorsD = eigenvectors + n * n;
        
        SharedSpace *Sharedcoef = new SharedSpace(3, num, 1, 1);
        double *coef = Sharedcoef->ObtainWriteEntireData();
        
		for (integer i = 0; i < num; i++)
		{
			dcopy_(&length, const_cast<double *> (L), &GLOBAL::IONE, Ltmp, &GLOBAL::IONE);
			/*Solve the linear system Li X = Lx, i.e., X = Li^{-1} Lx. The solution X is stored in LiiLx.
			Note that Li is a lower triangular matrix.
			Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
			dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, Ls + n * n * i, &N, Ltmp, &N, &info);
			if (info != 0)
			{
				std::cout << "The cholesky decompsotion in SPDMean::f failed with info:" << info << "!" << std::endl;
			}
			dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, Ltmp, &N, Ltmp, &N, &GLOBAL::DZERO, logLXL + n * n * i, &N);
            
             Matrix E(eigenvalues + n * i, n, 1), V(eigenvectors, n, n), VD(eigenvectorsD, n, n);
             Matrix MMt(logLXL + n * n * i, n, n);
             Matrix::EigenSymmetricM(GLOBAL::L, MMt, E, V);

             dcopy_(&length, eigenvectors, &GLOBAL::IONE, VD.matrix, &GLOBAL::IONE);
             for (integer j = 0; j < n; j++)
             {
                 a = log(eigenvalues[j + n * i]);
                 dscal_(&N, &a, eigenvectors + j * n, &GLOBAL::IONE);
             }
             
             
             Matrix MlogLXL(logLXL + n * n * i, n, n);
             Matrix::DGEMM(GLOBAL::DONE, V, false, VD, true, GLOBAL::DZERO, MlogLXL);

             
             coef[i] = 0.0;
             for(integer j = 0; j < N; j++) {coef[i] = coef[i] + log(eigenvalues[j + i * N]) * log(eigenvalues[j + i * N]);}
             coef[i] = std::sqrt(coef[i]);
             result = result + coef[i];
		}
        
        x->AddToTempData("logLXL", SharedlogLXL);
        x->AddToTempData("coef", Sharedcoef);

		delete[] Ltmp;
        delete[] eigenvalues;
        delete[] eigenvectors;
    

		return result;
	};
    
    

	void SPDMeanL1::RieGrad(Variable *x, Vector *gf) const
	{
		const SharedSpace *SharedlogLXL = x->ObtainReadTempData("logLXL");
		const double *logLXL = SharedlogLXL->ObtainReadData();
        const double *xM = x->ObtainReadData();
        
		double *gfVT = gf->ObtainWriteEntireData();
		for (integer i = 0; i < n * n; i++)
			gfVT[i] = 0;
		

		double *tmp = new double[n * n];
		integer N = n, info;
        
        
        const SharedSpace *Sharedcoef = x->ObtainReadTempData("coef");
        const double *coef  = Sharedcoef->ObtainReadData();
        
        double coefa, minval;
        
        
        minval = coef[0];
        integer idx = 0;
        for (integer i = 1; i < num; i++)
        {
            if(coef[i] < minval)
            {
                minval = coef[i];
                idx = i;
            }
        }
        
        
        
        if (minval > 1e-16)
        {
		    for (integer i = 0; i < num; i++)
		    {
			    /*tmp <-- log(Li^{-1} X Li^{-T}) Li^T */
			    dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, const_cast<double *> (logLXL + n * n * i), &N,
				    Ls + n * n * i, &N, &GLOBAL::DZERO, tmp, &N);

			    /*Solve the linear system Li^T X = tmp, i.e., X = Li^{-T} log(Li^{-1} X Li^{-T}) Li^T. The solution X is stored in tmp.
			    Note that Li is a lower triangular matrix.
			    Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
			    dtrtrs_(GLOBAL::L, GLOBAL::T, GLOBAL::N, &N, &N, Ls + n * n * i, &N, tmp, &N, &info);
			    if (info != 0)
			    {
				    std::cout << "The cholesky decompsotion in SPDMean::RieGrad failed with info:" << info << "!" << std::endl;
			    }
                coefa = 1.0/coef[i];
                dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &coefa, const_cast<double *> (xM), &N, tmp, &N, &GLOBAL::DONE, gfVT, &N);
		    }
		 }
		else
		 {
		    dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, const_cast<double *> (logLXL + n * n * idx), &N, 
		        Ls + n * n * idx, &N, &GLOBAL::DZERO, tmp, &N);
            
//          Li^{-T} * tmp -> tmp, i.e.,  Li^{-T} * logLXL * Li^T -> tmp
            dtrtrs_(GLOBAL::L, GLOBAL::T, GLOBAL::N, &N, &N, Ls + n * n * idx, &N, tmp, &N, &info);
            
            if (info != 0)
            {
                printf("The cholesky decompsotion in SPDMean::RieGrad failed with info:%d!\n", info);
            }
//          X * tmp -> eta
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *> (xM), &N, tmp, &N, &GLOBAL::DZERO, gfVT, &N);
		    
		 }
		
		
		delete[] tmp;
	};
    
    
    
    
    
//     double SPDMeanL1::dis(Variable *x) const
//     {
//          SPDManifold *Mani = dynamic_cast<SPDManifold *> (Domain);
//         // dis(x, y) = |log(Lx^{-1} y Lx^{-T})|_F
//         if (!x->TempDataExist("L"))
// 		{
// 			Mani->CholeskyRepresentation(x);
// 		}
// 		const SharedSpace *SharedLx = x->ObtainReadTempData("L");
// 		Variable *LElementx = SharedLx->GetSharedElement();
// 		const double *Lx = LElementx->ObtainReadData();
//          
//         integer N = n, info, length = n * n;
//         double *y = new double[n * n];
//         dcopy_(&length, Xtrue, &GLOBAL::IONE, y, &GLOBAL::IONE);
//          
//         dpotrf_(GLOBAL::L, &N, y, &N, &info);
//         for (integer j = 0; j < n; j++)
//         {
//             for (integer k = j + 1; k < n; k++)
//             {
//                 y[j + k * n] = 0; //manually make the upper triangle 0
//             }
//         }
// 
// 		double *LxLy = new double[n * n];
// //         
// 		dcopy_(&length, y, &GLOBAL::IONE, LxLy, &GLOBAL::IONE);
// // 		/*Solve the linear system Lx X = Ly, i.e., X = Lx^{-1}Ly. The solution X is stored in LxLy
// // 		Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
//         dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (Lx), &N, LxLy, &N, &info);
// //         
// //         /* compute temp = LxLy*(LxLy)^T  */
//         double *temp = new double[n * n];
//         dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, LxLy, &N, LxLy, &N, &GLOBAL::DZERO, temp, &N);
// //         
//         double *eigenvalues = new double[n + n * n];
// 	    double *eigenvectors = eigenvalues + n;
//         Matrix E(eigenvalues, n, 1), V(eigenvectors, n, n);
//         Matrix MMt(temp, n, n);
// 	    Matrix::EigenSymmetricM(GLOBAL::L, MMt, E, V);
//         
//         double result = 0.0;
//         for (integer i = 0; i < N; i++) 
//         {
//             result = result + log(eigenvalues[i])*log(eigenvalues[i]);
//         }
// 
//         result = std::sqrt(result);
// //         double result = 1.0;
//         
//         delete[] eigenvalues;   
//         delete[] temp;
//         delete[] y;
//         delete[] LxLy;
//         
//         return result;
//        
//     };
//     
//    
//     
//     
//     
//     double SPDMeanL1::distwo(Variable *x, Variable *y) const
//     {
//         SPDManifold *Mani = dynamic_cast<SPDManifold *> (Domain);
//         // dis(x, y) = |log(Lx^{-1} y Lx^{-T})|_F
//         if (!x->TempDataExist("L"))
// 		{
// 			Mani->CholeskyRepresentation(x);
// 		}
// 		const SharedSpace *SharedLx = x->ObtainReadTempData("L");
// 		Variable *LElementx = SharedLx->GetSharedElement();
// 		const double *Lx = LElementx->ObtainReadData();
//         
//         if (!y->TempDataExist("L"))
// 		{
// 			Mani->CholeskyRepresentation(y);
// 		}
// 		const SharedSpace *SharedLy = y->ObtainReadTempData("L");
// 		Variable *LElementy = SharedLy->GetSharedElement();
// 		const double *Ly = LElementy->ObtainReadData();
//         
//         integer N = n, info;
// 		double *LxLy = new double[N * N];
//         
//         integer length = n * n;
// 		dcopy_(&length, const_cast<double *> (Ly), &GLOBAL::IONE, LxLy, &GLOBAL::IONE);
// 		/*Solve the linear system Lx X = Ly, i.e., X = Lx^{-1}Ly. The solution X is stored in LxLy
// 		Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
//         dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (Lx), &N, LxLy, &N, &info);
//         
//         /* compute temp = LxLy*(LxLy)^T  */
//         double *temp = new double[N * N];
//         dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, LxLy, &N, LxLy, &N, &GLOBAL::DZERO, temp, &N);
//         
//         double *eigenvalues = new double[n + n * n];
// 	    double *eigenvectors = eigenvalues + n;
//         Matrix E(eigenvalues, n, 1), V(eigenvectors, n, n);
//         Matrix MMt(temp, n, n);
// 	    Matrix::EigenSymmetricM(GLOBAL::L, MMt, E, V);
//         
//         double result = 0.0;
//         for (integer i = 0; i < N; i++) 
//         {
//             result = result + log(eigenvalues[i])*log(eigenvalues[i]);
//         }
// 
//         result = sqrt(result);
//         return result;
//         
//         delete[] eigenvalues;   
//         delete[] temp;
//     }
//     
    

	void SPDMeanL1::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const
	{
		std::cout << "warning: SPDMean::RieHessianEta has not been implemented!" << std::endl;
		etax->CopyTo(xix);
	};
}; /*end of ROPTLIB namespace*/
