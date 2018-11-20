
#include "SPDMeanStein.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDMeanStein::SPDMeanStein(double *inLs, double *inXtrue, integer inn, integer innum)
	{
		Ls = inLs;
        Xtrue = inXtrue;
		n = inn;
		num = innum;
	};

	SPDMeanStein::~SPDMeanStein(void)
	{
	};

	double SPDMeanStein::f(Variable *x) const
	{
        SPDManifold *Mani = dynamic_cast<SPDManifold *> (Domain);
        if (!x->TempDataExist("L"))
        {
            Mani->CholeskyRepresentation(x);
        }
        const SharedSpace *SharedL = x->ObtainReadTempData("L");
        Variable *LElement = SharedL->GetSharedElement();
        const double *L = LElement->ObtainReadData();
        
        
        const double *xxM = x->ObtainReadData();
        
        integer N = n, length = n * n, Nnum = n * n * num;
        double onehalf = 0.5;
        
        double *I = new double[n * n]; // make identity matrix
        for (integer i = 0; i < length; i++) I[i] = 0.0;
        for (integer i = 0; i < n; i++) I[i * n + i] = 1.0;
        
        
//         for (integer i = 0; i < num; i++)
//         {
//         printf("i: %d", i);
//         ForDebug::Print("Ls:", Ls + n * n * i, n, n);
//         }
//         
        
        // Used to store Ai*x
        SharedSpace *SharedLplusX = new SharedSpace(3, n, n, num);
        double *LplusX = SharedLplusX->ObtainWriteEntireData();
        
        double *eigenvalues = new double[n + n * n];
        double *eigenvectors = eigenvalues + n;
        double *tmp = new double[n * n];
        for (integer i = 0; i < length; i++) tmp[i] = 0.0;
        
       
        dcopy_(&Nnum, Ls, &GLOBAL::IONE, LplusX, &GLOBAL::IONE);
        double *coef = new double[num]; // used to store trace(log(Ai + X)/2)
        double *coef1 = new double[num];
        double *coef2 = new double[num]; // used to store trace(log(Ai + X)/2)
        double *tmp1 = new double[n * n];
        
        
//        for(integer i = 0; i < num; i++)
//            ForDebug::Print("Ls:", Ls + n * n * i, n, n);
        
        
//         ForDebug::Print("X:", const_cast<double *>(xxM), n, n);
        
        
        for (integer i = 0; i < num; i++)
        {
            // compute 0.5 * (Ai + x)
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &onehalf, const_cast<double *>(xxM), &N, I, &N, &onehalf, LplusX + n * n * i, &N);
            

//            printf("i: %d", i);
//            ForDebug::Print("LplusX:", LplusX + n * n * i, n, n);
  

            // eigenvalue decomposition of LplusX
            Matrix E(eigenvalues, n, 1), V(eigenvectors, n, n);
            Matrix MMt(LplusX + n * n * i, n, n);
            Matrix::EigenSymmetricM(GLOBAL::L, MMt, E, V);
            
            coef1[i] = 0.0;
            for (integer j = 0; j < n; j++)
            {
//                printf("eigenvalues: %.4e\n", eigenvalues[j]);
                coef1[i] = coef1[i] + log(eigenvalues[j]);
            }
            
            
            // bug here !!!!!!!!!!!
            
            
            // compute L^T * Ai * L
            
            // compute L^T * Ai -> tmp
            dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *>(L), &N, Ls + n * n * i, &N, &GLOBAL::DZERO, tmp, &N);
            
            
            // compute tmp * L -> tmp
            dcopy_(&length, tmp, &GLOBAL::IONE, tmp1, &GLOBAL::IONE);
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp1, &N, const_cast<double *>(L), &N, &GLOBAL::DZERO, tmp, &N);
            
//             printf("i: %d", i);
//             ForDebug::Print("Ax:", tmp, n, n);
            
            
            // eigenvalue decomposition of LplusX
            Matrix E2(eigenvalues, n, 1), V2(eigenvectors, n, n);
            Matrix MMt2(tmp, n, n);
            Matrix::EigenSymmetricM(GLOBAL::L, MMt2, E2, V2);
            
            
            coef2[i] = 0.0;
            for (integer j = 0; j < n; j++)
            {
//                 printf("eigenvalues: %.4e\n", eigenvalues[j]);
                coef2[i] = coef2[i] + 0.5 * log(eigenvalues[j]);
            
            }
            

            for (integer i = 0; i < num; i++) coef[i] = coef1[i] - coef2[i];
            
            

        }
        x->AddToTempData("LplusX", SharedLplusX);
        


        double result = 0.0;
        for (integer i = 0; i < num; i++)
        {
            result = result + coef[i];
        }
        
//         for (integer i = 0; i < num; i++) printf("coef1: %.4e\n", coef1[i]);
//         for (integer i = 0; i < num; i++) printf("coef2: %.4e\n", coef2[i]);
//         for (integer i = 0; i < num; i++) printf("coef: %.4e\n", coef[i]);
        

		result = result / (2.0 * num);
        
        
        delete[] eigenvalues;
        delete[] I;
        delete[] tmp;
        delete[] tmp1;
        delete[] coef1;
        delete[] coef2;
        delete[] coef;
        
        
		return result;
	};
    
    

    void SPDMeanStein::RieGrad(Variable *x, Vector *gf) const
    {
        const SharedSpace *SharedLplusX = x->ObtainReadTempData("LplusX");
        const double *LplusX = SharedLplusX->ObtainReadData();
        
        const double *xM = x->ObtainReadData();
        
        double *gfVT = gf->ObtainWriteEntireData();
        
        for (integer i = 0; i < n * n; i++)
            gfVT[i] = 0.0;
        
        
        double *tmp = new double[n * n];
        double *tmp1 = new double[n * n];
        double *LLplusX = new double[n * n * num];
        
        SharedSpace *SharediLplusX = new SharedSpace(3, n, n, num);
        double *iLplusX = SharediLplusX->ObtainWriteEntireData();
        
        integer length = n * n, N = n, Nnum = n * n * num, info;
        double negone = -1.0;
        
        double *I = new double[n * n]; // make identity matrix
        for (integer i = 0; i < length; i++) I[i] = 0.0;
        for (integer i = 0; i < n; i++) I[i * n + i] = 1.0;
        
        
        dcopy_(&Nnum, const_cast<double *>(LplusX), &GLOBAL::IONE, LLplusX, &GLOBAL::IONE);
        
        
        for (integer i = 0; i < num; i++)
        {
            // Cholesky factorization of LplusX -> iLplusX
            dpotrf_(GLOBAL::L, &N, LLplusX + n * n * i, &N, &info);
            for (integer j = 0; j < n; j++)
                for (integer k = j + 1; k < n; k++)
                    LLplusX[j + k * n + i * n * n] = 0;
            
            // Compute inv(LplusX) -> iLplusX
            dcopy_(&length, I, &GLOBAL::IONE, iLplusX + n * n * i, &GLOBAL::IONE);
            dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, LLplusX + n * n * i, &N, iLplusX + n * n * i, &N, &info);
            
            // compute inv(LplusX) = inv(LplusX)^T * inv(LplusX)
            dcopy_(&length, iLplusX + n * n * i, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, tmp, &N, &GLOBAL::DZERO, iLplusX + n * n * i, &N);
            
            // compute X * inv(LplusX) -> tmp
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *>(xM), &N, iLplusX + n * n * i, &N, &GLOBAL::DZERO, tmp, &N);
            
            // compute iLplusX * X - X -> tmp1
            dcopy_(&length, const_cast<double *>(xM), &GLOBAL::IONE, tmp1, &GLOBAL::IONE);
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, const_cast<double *>(xM), &N, &negone, tmp1, &N);
            
            
            for(integer j = 0; j < n * n; j++) gfVT[j] = gfVT[j] + tmp1[j];
            
        }
        
        
        double scalar = 1.0/(4.0 * num);
        dscal_(&length, &scalar, gfVT, &GLOBAL::IONE);
        
        x -> AddToTempData("iLplusX", SharediLplusX);
        
        delete[] I;
        delete[] tmp;
        delete[] tmp1;
        delete[] LLplusX;
    };
    
    
    void SPDMeanStein::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const
    {
        
        if (!x -> TempDataExist("iLplusX"))
        {
            const double *xM = x->ObtainReadData();
            
            const SharedSpace *Temp = x->ObtainReadTempData("LplusX");
            const double *LplusX = Temp->ObtainReadData();
            
            SharedSpace *iLplusX = new SharedSpace(3, n, n, num);
            double *iLplusXptr = iLplusX->ObtainWriteEntireData();
            
            integer N = n, length = N * N, numlength = N * N * num;
            integer info;
            
            double *I = new double[n * n]; // make identity matrix
            for (integer i = 0; i < length; i++) I[i] = 0.0;
            for (integer i = 0; i < n; i++) I[i * n + i] = 1.0;
            
            double *LLplusX = new double[n * n * num];
            double *tmp = new double[n * n];
            
            dcopy_(&numlength, const_cast<double *>(LplusX), &GLOBAL::IONE, LLplusX, &GLOBAL::IONE);
            
            
            for (integer i = 0; i < num; i++)
            {
                // Cholesky factorization of LplusX -> iLplusX
                dpotrf_(GLOBAL::L, &N, LLplusX + n * n * i, &N, &info);
                for (integer j = 0; j < n; j++)
                    for (integer k = j + 1; k < n; k++)
                        LLplusX[j + k * n + n * n * i] = 0;
                
                // Compute inv(LplusX) -> iLplusX
                dcopy_(&length, I, &GLOBAL::IONE, iLplusXptr + n * n * i, &GLOBAL::IONE);
                dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, LLplusX + n * n * i, &N, iLplusXptr + n * n * i, &N, &info);
                //            ForDebug::Print("inv(LplusX):", iLplusX + n * n * i, n, n);
                
                // compute inv(LplusX) = inv(LplusX)^T * inv(LplusX)
                dcopy_(&length, iLplusXptr + n * n * i, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
                dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, tmp, &N, &GLOBAL::DZERO, iLplusXptr + n * n * i, &N);
            }
            x->AddToTempData("iLplusX", iLplusX);
            
            delete[] I;
            delete[] LLplusX;
            delete[] tmp;
        }
        
        integer N = n, length = N * N, numlength = num * length;
        
        const SharedSpace *SharediLplusX = x->ObtainReadTempData("iLplusX");
        const double *iLplusX = SharediLplusX->ObtainReadData();
        const double *xM = x->ObtainReadData();
        const double *etaxTV = etax->ObtainReadData();
        double *xixTV = xix->ObtainWriteEntireData();
        
        
        
        
        if (xixTV == etaxTV)
        {
            std::cout << "Error in RieHessianEta!" << std::endl;
            exit(0);
        }
        
        double *tmp1 = new double[n * n];
        double *tmp2 = new double[n * n];
        double *tmp3 = new double[n * n];
        double *rightX = new double[n * n];
        double *leftX = new double[n * n];
        
        
        for (integer i = 0; i < n * n; i++)
            xixTV[i] = 0.0;
        
        for (integer i = 0; i < num; i++)
        {
            // iL * X -> rightX
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *>(iLplusX + n * n * i), &N, const_cast<double *>(xM), &N, &GLOBAL::DZERO, rightX, &N);
            
            // X * iL -> leftX
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *>(xM), &N, const_cast<double *>(iLplusX + n * n * i), &N, &GLOBAL::DZERO, leftX, &N);
            
            // xi * iL * X -> tmp1
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *>(etaxTV), &N, rightX, &N, &GLOBAL::DZERO, tmp1, &N);
            
            // X * iL * xi -> tmp2
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, leftX, &N, const_cast<double *>(etaxTV), &N, &GLOBAL::DZERO, tmp2, &N);
            
            // tmp3
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp2, &N, rightX, &N, &GLOBAL::DZERO, tmp3, &N);
            
            for(integer j = 0; j < n * n; j++) xixTV[j] = xixTV[j] + tmp1[j] + tmp2[j] - tmp3[j];
            
        }
        
        
        double scalar = 1.0/(8.0 * num);
        dscal_(&length, &scalar, xixTV, &GLOBAL::IONE);
        
        delete[] tmp1;
        delete[] tmp2;
        delete[] tmp3;
        delete[] rightX;
        delete[] leftX;
    };
    
    
    
    double SPDMeanStein::dis(Variable *x) const
    {
         SPDManifold *Mani = dynamic_cast<SPDManifold *> (Domain);
        // dis(x, y) = |log(Lx^{-1} y Lx^{-T})|_F
        if (!x->TempDataExist("L"))
		{
			Mani->CholeskyRepresentation(x);
		}
		const SharedSpace *SharedLx = x->ObtainReadTempData("L");
		Variable *LElementx = SharedLx->GetSharedElement();
		const double *Lx = LElementx->ObtainReadData();
         
        integer N = n, info, length = n * n;
        double *y = new double[n * n];
        dcopy_(&length, Xtrue, &GLOBAL::IONE, y, &GLOBAL::IONE);
         
        dpotrf_(GLOBAL::L, &N, y, &N, &info);
        for (integer j = 0; j < n; j++)
        {
            for (integer k = j + 1; k < n; k++)
            {
                y[j + k * n] = 0; //manually make the upper triangle 0
            }
        }

		double *LxLy = new double[n * n];
//         
		dcopy_(&length, y, &GLOBAL::IONE, LxLy, &GLOBAL::IONE);
// 		/*Solve the linear system Lx X = Ly, i.e., X = Lx^{-1}Ly. The solution X is stored in LxLy
// 		Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (Lx), &N, LxLy, &N, &info);
//         
//         /* compute temp = LxLy*(LxLy)^T  */
        double *temp = new double[n * n];
        dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, LxLy, &N, LxLy, &N, &GLOBAL::DZERO, temp, &N);
//         
        double *eigenvalues = new double[n + n * n];
	    double *eigenvectors = eigenvalues + n;
        Matrix E(eigenvalues, n, 1), V(eigenvectors, n, n);
        Matrix MMt(temp, n, n);
	    Matrix::EigenSymmetricM(GLOBAL::L, MMt, E, V);
        
        double result = 0.0;
        for (integer i = 0; i < N; i++) 
        {
            result = result + log(eigenvalues[i])*log(eigenvalues[i]);
        }

        result = std::sqrt(result);
//         double result = 1.0;
        
        delete[] eigenvalues;   
        delete[] temp;
        delete[] y;
        delete[] LxLy;
        
        return result;
       
    };
    
    
    
    
    double SPDMeanStein::distwo(Variable *x, Variable *y) const
    {
        SPDManifold *Mani = dynamic_cast<SPDManifold *> (Domain);
        // dis(x, y) = |log(Lx^{-1} y Lx^{-T})|_F
        if (!x->TempDataExist("L"))
		{
			Mani->CholeskyRepresentation(x);
		}
		const SharedSpace *SharedLx = x->ObtainReadTempData("L");
		Variable *LElementx = SharedLx->GetSharedElement();
		const double *Lx = LElementx->ObtainReadData();
        
        if (!y->TempDataExist("L"))
		{
			Mani->CholeskyRepresentation(y);
		}
		const SharedSpace *SharedLy = y->ObtainReadTempData("L");
		Variable *LElementy = SharedLy->GetSharedElement();
		const double *Ly = LElementy->ObtainReadData();
        
        integer N = n, info;
		double *LxLy = new double[N * N];
        
        integer length = n * n;
		dcopy_(&length, const_cast<double *> (Ly), &GLOBAL::IONE, LxLy, &GLOBAL::IONE);
		/*Solve the linear system Lx X = Ly, i.e., X = Lx^{-1}Ly. The solution X is stored in LxLy
		Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (Lx), &N, LxLy, &N, &info);
        
        /* compute temp = LxLy*(LxLy)^T  */
        double *temp = new double[N * N];
        dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, LxLy, &N, LxLy, &N, &GLOBAL::DZERO, temp, &N);
        
        double *eigenvalues = new double[n + n * n];
	    double *eigenvectors = eigenvalues + n;
        Matrix E(eigenvalues, n, 1), V(eigenvectors, n, n);
        Matrix MMt(temp, n, n);
	    Matrix::EigenSymmetricM(GLOBAL::L, MMt, E, V);
        
        double result = 0.0;
        for (integer i = 0; i < N; i++) 
        {
            result = result + log(eigenvalues[i])*log(eigenvalues[i]);
        }

        result = sqrt(result);
        
        delete[] eigenvalues;   
        delete[] temp;
        return result;
        
    }
    
    
    
    
}; /*end of ROPTLIB namespace*/
