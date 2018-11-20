
#include "SPDMeanSymmLDOneParamLinfty.h"

/*Define the namespace*/
namespace ROPTLIB{
    
    SPDMeanSymmLDOneParamLinfty::SPDMeanSymmLDOneParamLinfty(double *inLs, double *inLDAs, integer inn, integer innum, double inalpha)
    {
        Ls = inLs;
        LDAs = inLDAs;
        n = inn;
        num = innum;
        alpha = inalpha;
    };
    
    SPDMeanSymmLDOneParamLinfty::~SPDMeanSymmLDOneParamLinfty(void)
    {
    };
    
    double SPDMeanSymmLDOneParamLinfty::f(Variable *x) const
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
        
        // make identity matrix
        double *I = new double[n * n];
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
        
        
        
        
        
        // XplusCL is used to store a1 * x + a2 * Ai
        SharedSpace *SharedXplusCL = new SharedSpace(3, n, n, num);
        double *XplusCL = SharedXplusCL->ObtainWriteEntireData();
        double *coef2 = new double[num];
        dcopy_(&Nnum, Ls, &GLOBAL::IONE, XplusCL, &GLOBAL::IONE);
        
        for (integer i = 0; i < num; i++)
        {
            // compute a2 * Ai + a1 * x
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &a1, const_cast<double *>(xxM), &N, I, &N, &a2, XplusCL + n * n * i, &N);
            
            // Cholesky factorization of XplusCL
            dpotrf_(GLOBAL::L, &N, XplusCL + n * n * i, &N, &info);
            for (integer j = 0; j < n; j++)
                for (integer k = j + 1; k < n; k++)
                    XplusCL[j + k * n + i * n * n] = 0;
            
            // compute the determinant
            coef2[i] = 1.0;
            for (integer j = 0; j < n; j++)
            {
                coef2[i] = coef2[i] * XplusCL[j * n + j + i * n * n];
            }
            coef2[i] = log(coef2[i] * coef2[i]);
        }
        x->AddToTempData("XplusCL", SharedXplusCL);
        
        
        
        // compute the determinant of X: product of diagonal elements
        double detX = 1.0;
        for (integer j = 0; j < n; j++)
        {
            detX = detX * const_cast<double *>(L)[j * n + j];
        }
        detX = log(detX * detX);
        
        
        SharedSpace *Sharedmaxval = new SharedSpace(1, 2); // maxval is used to store the index and the max value
        double *maxval = Sharedmaxval->ObtainWriteEntireData();
        
        
        for (integer i = 0; i < num; i++)
        {
            coef1[i] = coef1[i] + coef2[i] - LDAs[i] - detX;
        }
        
        
        double result = coef1[0];
        maxval[1] = 0;
        for (integer i = 1; i < num; i++)
        {
            if(coef1[i] > result)
            {
                result = coef1[i];
                maxval[1] = static_cast<double>(i);
            }
        }
        maxval[0] = result;
        
        result = (result * 2.0)/(1 - alpha * alpha);
        x->AddToTempData("maxval", Sharedmaxval);
        
        
        delete[] I;
        delete[] coef1;
        delete[] coef2;
        
        return result;
        
    };
    
    
    
    void SPDMeanSymmLDOneParamLinfty::RieGrad(Variable *x, Vector *gf) const
    {
        const SharedSpace *Sharedmaxval = x->ObtainReadTempData("maxval");
        const double *maxval  = Sharedmaxval->ObtainReadData();
        
        const SharedSpace *SharedCLplusX = x->ObtainReadTempData("CLplusX");
        const double *CLplusX = SharedCLplusX->ObtainReadData();
        
        const SharedSpace *SharedXplusCL = x->ObtainReadTempData("XplusCL");
        const double *XplusCL = SharedXplusCL->ObtainReadData();
        
        const double *xM = x->ObtainReadData();
        double *gfVT = gf->ObtainWriteEntireData();
        for (integer i = 0; i < n * n; i++) gfVT[i] = 0.0;
        
        double *tmp = new double[n * n];
        double *iLplusX = new double[n * n];
        double *iXplusL = new double[n * n];
        
        
        integer length = n * n, N = n, Nnum = n * n * num, info;
        double negone = -1.0, a1 = (1.0 - alpha)/2.0, a2 = (1.0 + alpha)/2.0;
        integer idx = static_cast<integer>(maxval[1]);

        
        double *I = new double[n * n]; // make identity matrix
        for (integer i = 0; i < length; i++) I[i] = 0.0;
        for (integer i = 0; i < n; i++) I[i * n + i] = 1.0;
        
        
               
        // Compute inv(LplusX) -> iLplusX
        dcopy_(&length, const_cast<double *>(CLplusX + n * n * idx), &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dcopy_(&length, I, &GLOBAL::IONE, iLplusX, &GLOBAL::IONE);
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, tmp, &N, iLplusX, &N, &info);
        
        // compute inv(LplusX) = inv(LplusX)^T * inv(LplusX) -> iLplusX
        dcopy_(&length, iLplusX, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, tmp, &N, &GLOBAL::DZERO, iLplusX, &N);
        
        
        // Compute inv(XplusL) -> iXplusL
        dcopy_(&length, const_cast<double *>(XplusCL + n * n * idx), &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dcopy_(&length, I, &GLOBAL::IONE, iXplusL, &GLOBAL::IONE);
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, tmp, &N, iXplusL, &N, &info);
        
        // compute inv(LplusX) = inv(LplusX)^T * inv(LplusX)
        dcopy_(&length, iXplusL, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, tmp, &N, &GLOBAL::DZERO, iXplusL, &N);
        
        
        // a2 * iLplusX + a1 * iXplusL -> gfVT
        dcopy_(&length, iXplusL, &GLOBAL::IONE, gfVT, &GLOBAL::IONE);
        dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &a2, iLplusX, &N, I, &N, &a1, gfVT, &N);
        
        
        // compute X * gfVT -> tmp
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *>(xM), &N, gfVT, &N, &GLOBAL::DZERO, tmp, &N);
        
        // compute tmp * X -> gfVT
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, const_cast<double *>(xM), &N, &GLOBAL::DZERO, gfVT, &N);
        
        // compute gfVT - X -> tmp
        dcopy_(&length, gfVT, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dcopy_(&length, const_cast<double *>(xM), &GLOBAL::IONE, gfVT, &GLOBAL::IONE);
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, I, &N, &negone, gfVT, &N);
        
        double scalar = 2.0/(1.0 - alpha);
        dscal_(&length, &scalar, gfVT, &GLOBAL::IONE);
        
        
        delete[] I;
        delete[] tmp;
        delete[] iLplusX;
        delete[] iXplusL;
        
        
    };
    
    
    
    void SPDMeanSymmLDOneParamLinfty::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const
    {
        std::cout << "warning: SPDMeanSymmLDOneParam::RieHessianEta has not been implemented!" << std::endl;
        etax->CopyTo(xix);
    };
}; /*end of ROPTLIB namespace*/
