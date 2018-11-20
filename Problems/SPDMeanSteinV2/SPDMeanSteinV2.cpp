
#include "SPDMeanSteinV2.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDMeanSteinV2::SPDMeanSteinV2(double *inLs, integer inn, integer innum)
	{
		Ls = inLs;
		n = inn;
		num = innum;
	};

	SPDMeanSteinV2::~SPDMeanSteinV2(void)
	{
	};


    double SPDMeanSteinV2::f(Variable *x) const
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
        
        integer N = n, length = n * n, Nnum = n * n * num, info;
        double onehalf = 0.5;
        
        double *I = new double[n * n]; // make identity matrix
        for (integer i = 0; i < length; i++) I[i] = 0.0;
        for (integer i = 0; i < n; i++) I[i * n + i] = 1.0;
        
        
        // CLplusX is used to store Ai + x        
        SharedSpace *SharedCLplusX = new SharedSpace(3, n, n, num);
        double *CLplusX = SharedCLplusX->ObtainWriteEntireData();
        
        double *tmp = new double[n * n];
        double *tmp1 = new double[n * n];
        double *coef1 = new double[num];
        double *coef2 = new double[num]; // used to store trace(log(Ai + X)/2)
        
       
        dcopy_(&Nnum, Ls, &GLOBAL::IONE, CLplusX, &GLOBAL::IONE);
        
        for (integer i = 0; i < num; i++)
        {
            // compute 0.5 * (Ai + x)
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &onehalf, const_cast<double *>(xxM), &N, I, &N, &onehalf, CLplusX + n * n * i, &N);
            
            // Cholesky factorization of CLplusX
            dpotrf_(GLOBAL::L, &N, CLplusX + n * n * i, &N, &info);
            for (integer j = 0; j < n; j++)
                for (integer k = j + 1; k < n; k++)
                    CLplusX[j + k * n + i * n * n] = 0;
            
            // compute the determinant 
            coef1[i] = 1.0;
            for (integer j = 0; j < n; j++)
            {
                coef1[i] = coef1[i] * CLplusX[j * n + j + i * n * n];
            }
            coef1[i] = log(coef1[i] * coef1[i]);            
            
            
            // compute L^T * Ai * L -> tmp
            // compute L^T * Ai -> tmp
            dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *>(L), &N, Ls + n * n * i, &N, &GLOBAL::DZERO, tmp, &N);
            
            // compute tmp * L -> tmp
            dcopy_(&length, tmp, &GLOBAL::IONE, tmp1, &GLOBAL::IONE);
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp1, &N, const_cast<double *>(L), &N, &GLOBAL::DZERO, tmp, &N);
            
            // Cholesky factorization of L^T * Ai * L
            dpotrf_(GLOBAL::L, &N, tmp, &N, &info);
            
            // compute determinant: product of diagonal elements
            coef2[i] = 1.0;
            for (integer j = 0; j < n; j++)
            {
                coef2[i] = coef2[i] * tmp[j * n + j];
            }
            coef2[i] = 0.5 * log(coef2[i] * coef2[i]);
            
        }
        x->AddToTempData("CLplusX", SharedCLplusX);

        
        double result = 0.0;
        for (integer i = 0; i < num; i++)
        {
            result = result + coef1[i] - coef2[i];
        }

        result = result / (2.0 * num);
        
        
        delete[] I;
        delete[] tmp;
        delete[] tmp1;
        delete[] coef1;
        delete[] coef2;

        return result;
    };
    
    
    
    void SPDMeanSteinV2::RieGrad(Variable *x, Vector *gf) const
    {        
        const SharedSpace *SharedCLplusX = x->ObtainReadTempData("CLplusX");
        const double *CLplusX = SharedCLplusX->ObtainReadData();
        
        SharedSpace *SharediLplusX = new SharedSpace(3, n, n, num); // used to store the inverse of L
        double *iLplusX = SharediLplusX->ObtainWriteEntireData();
        
        
        const double *xM = x->ObtainReadData();

        double *gfVT = gf->ObtainWriteEntireData();
        
        for (integer i = 0; i < n * n; i++)
            gfVT[i] = 0.0;
        
        
        double *tmp = new double[n * n];
        double *tmp1 = new double[n * n];
        double *LLplusX = new double[n * n * num];
        
        integer length = n * n, N = n, Nnum = n * n * num, info;
        double negone = -1.0;
        double *I = new double[n * n]; // make identity matrix
        for (integer i = 0; i < length; i++) I[i] = 0.0;
        for (integer i = 0; i < n; i++) I[i * n + i] = 1.0;

        
        dcopy_(&Nnum, const_cast<double *>(CLplusX), &GLOBAL::IONE, LLplusX, &GLOBAL::IONE);
        
        for (integer i = 0; i < num; i++)
        {
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
    
    
    void SPDMeanSteinV2::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const
    {
        
        if (!x -> TempDataExist("iLplusX"))
        {
            const double *xM = x->ObtainReadData();
            
            const SharedSpace *Temp = x->ObtainReadTempData("CLplusX");
            const double *LplusX = Temp->ObtainReadData();
            
            SharedSpace *iLplusX = new SharedSpace(3, n, n, num);
            double *iLplusXptr = iLplusX->ObtainWriteEntireData();
            
            integer N = n, length = N * N, numlength = N * N * num, info;
            
            double *I = new double[n * n]; // make identity matrix
            for (integer i = 0; i < length; i++) I[i] = 0.0;
            for (integer i = 0; i < n; i++) I[i * n + i] = 1.0;
            
            double *LLplusX = new double[n * n * num];
            double *tmp = new double[n * n];
            
            dcopy_(&numlength, const_cast<double *>(LplusX), &GLOBAL::IONE, LLplusX, &GLOBAL::IONE);
            
            
            for (integer i = 0; i < num; i++)
            {
//                 // Cholesky factorization of LplusX -> iLplusX
//                 dpotrf_(GLOBAL::L, &N, LLplusX + n * n * i, &N, &info);
//                 for (integer j = 0; j < n; j++)
//                     for (integer k = j + 1; k < n; k++)
//                         LLplusX[j + k * n + n * n * i] = 0;
                
                // Compute inv(LplusX) -> iLplusX
                dcopy_(&length, I, &GLOBAL::IONE, iLplusXptr + n * n * i, &GLOBAL::IONE);
                dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, LLplusX + n * n * i, &N, iLplusXptr + n * n * i, &N, &info);
                
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
  
}; /*end of ROPTLIB namespace*/
