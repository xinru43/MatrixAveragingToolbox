
#include "SPDMeanLDOneParam.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDMeanLDOneParam::SPDMeanLDOneParam(double *inLs, integer inn, integer innum, double inalpha)
	{
		Ls = inLs;
		n = inn;
		num = innum;
        alpha = inalpha;
	};

	SPDMeanLDOneParam::~SPDMeanLDOneParam(void)
	{
	};


    double SPDMeanLDOneParam::f(Variable *x) const
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
        double a1 = (1.0 - alpha)/2.0, a2 = (1.0 + alpha)/2.0;
        
        double *I = new double[n * n]; // make identity matrix
        for (integer i = 0; i < length; i++) I[i] = 0.0;
        for (integer i = 0; i < n; i++) I[i * n + i] = 1.0;
        
        
        // CLplusX is used to store a1 * Ai + a2 * x        
        SharedSpace *SharedCLplusX = new SharedSpace(3, n, n, num);
        double *CLplusX = SharedCLplusX->ObtainWriteEntireData();
        
        double *coef1 = new double[num];
        
       
        dcopy_(&Nnum, Ls, &GLOBAL::IONE, CLplusX, &GLOBAL::IONE);
        
        for (integer i = 0; i < num; i++)
        {
            // compute a1 * Ai + a2 * x
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &a2, const_cast<double *>(xxM), &N, I, &N, &a1, CLplusX + n * n * i, &N);
            
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
        }
        x->AddToTempData("CLplusX", SharedCLplusX);
        
        
        // compute the determinant of X: product of diagonal elements
        double detX = 1.0;
        for (integer j = 0; j < n; j++)
        {
            detX = detX * const_cast<double *>(L)[j * n + j];
        }
        detX = log(detX * detX);
        detX = a2 * num * detX;

        
        double result = 0.0;
        for (integer i = 0; i < num; i++)
        {
            result = result + coef1[i];
        }
        result = result - detX;
        result = (result * 4.0)/(1 - alpha * alpha);
        
        
        delete[] I;
        delete[] coef1;

        return result;
    };
    
    
    
    void SPDMeanLDOneParam::RieGrad(Variable *x, Vector *gf) const
    {        
        const SharedSpace *SharedCLplusX = x->ObtainReadTempData("CLplusX");
        const double *CLplusX = SharedCLplusX->ObtainReadData();
        
        SharedSpace *SharediLplusX = new SharedSpace(3, n, n, num); // used to store the inverse of L
        double *iLplusX = SharediLplusX->ObtainWriteEntireData();
        
        
        const double *xM = x->ObtainReadData();
        double *gfVT = gf->ObtainWriteEntireData();
        
        for (integer i = 0; i < n * n; i++) gfVT[i] = 0.0;
        
        
        double *tmp = new double[n * n];
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
            
            // add each iLplusX to gfVT
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, iLplusX + n * n * i, &N, I, &N, &GLOBAL::DONE, gfVT, &N);
        }
        
        
         // compute X * gfVT -> tmp
         dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *>(xM), &N, gfVT, &N, &GLOBAL::DZERO, tmp, &N);
            
        // compute tmp * X -> gfVT
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, const_cast<double *>(xM), &N, &GLOBAL::DZERO, gfVT, &N);

        // compute gfVT - num * X -> tmp
        dcopy_(&length, gfVT, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dcopy_(&length, const_cast<double *>(xM), &GLOBAL::IONE, gfVT, &GLOBAL::IONE);
        double negnum = -1.0 * num;
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, I, &N, &negnum, gfVT, &N);

        
        double scalar = 2.0/(1.0 - alpha);
        dscal_(&length, &scalar, gfVT, &GLOBAL::IONE);
        
        x -> AddToTempData("iLplusX", SharediLplusX);
        
        delete[] I;
        delete[] tmp;
        delete[] LLplusX;
    };
    
    
    
    
    
    
    void SPDMeanLDOneParam::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const
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
            
            for(integer j = 0; j < n * n; j++) xixTV[j] = xixTV[j] + tmp1[j] + tmp2[j] - (1.0 + alpha)*tmp3[j];
            
        }
        
        
        double scalar = 1.0/(1.0 - alpha);
        dscal_(&length, &scalar, xixTV, &GLOBAL::IONE);
        
        delete[] tmp1;
        delete[] tmp2;
        delete[] tmp3;
        delete[] rightX;
        delete[] leftX;
    };
  
}; /*end of ROPTLIB namespace*/
