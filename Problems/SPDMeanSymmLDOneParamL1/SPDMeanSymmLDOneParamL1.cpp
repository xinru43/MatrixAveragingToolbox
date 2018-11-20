
#include "SPDMeanSymmLDOneParamL1.h"

/*Define the namespace*/
namespace ROPTLIB{
    
    SPDMeanSymmLDOneParamL1::SPDMeanSymmLDOneParamL1(double *inLs, double *inLDAs, integer inn, integer innum, double inalpha)
    {
        Ls = inLs;
        LDAs = inLDAs;
        n = inn;
        num = innum;
        alpha = inalpha;
    };
    
    SPDMeanSymmLDOneParamL1::~SPDMeanSymmLDOneParamL1(void)
    {
    };
    
    
    double SPDMeanSymmLDOneParamL1::f(Variable *x) const
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
        
        // XplusCL is used to store a1 * x + a2 * Ai
        SharedSpace *SharedXplusCL = new SharedSpace(3, n, n, num);
        
        
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
        
//         // compute a1 * log(det A_i) -> lda
//         SharedSpace *Sharedlda = new SharedSpace(3, num, 1, 1);
//         double *lda = Sharedlda->ObtainWriteEntireData();
//         double *CLs = new double[n * n * num];
//         
//         dcopy_(&Nnum, Ls, &GLOBAL::IONE, CLs, &GLOBAL::IONE);
//         for (integer i = 0; i < num; i++)
//         {
//             // Cholesky factorization of Ls -> CLs
//             dpotrf_(GLOBAL::L, &N, CLs + n * n * i, &N, &info);
//             for (integer j = 0; j < n; j++)
//                 for (integer k = j + 1; k < n; k++)
//                     CLs[j + k * n + i * n * n] = 0;
//             
//             // compute the determinant
//             lda[i] = 1.0;
//             for (integer j = 0; j < n; j++)
//             {
//                 lda[i] = lda[i] * CLs[j * n + j + i * n * n];
//             }
//             lda[i] = log(lda[i] * lda[i]);
//             
//         }
//         x->AddToTempData("lda", Sharedlda);
        
        

        // compute the determinant of X: product of diagonal elements
        double detX = 1.0;
        for (integer j = 0; j < n; j++)
        {
            detX = detX * const_cast<double *>(L)[j * n + j];
        }
        detX = log(detX * detX);

        
        
        SharedSpace *Sharedcoef3 = new SharedSpace(3, num, 1, 1);
        double *coef3 = Sharedcoef3->ObtainWriteEntireData();
        
        double result = 0.0;
        for (integer i = 0; i < num; i++)
        {
            coef3[i] = coef1[i] + coef2[i] - detX - LDAs[i];
            coef3[i] = sqrt(coef3[i]);
            result = result + coef3[i];
        }
        x->AddToTempData("coef3", Sharedcoef3);
               
        delete[] I;
        delete[] coef1;
        delete[] coef2;
        
        return result;
    };
    
    
    
    
    
    
    void SPDMeanSymmLDOneParamL1::RieGrad(Variable *x, Vector *gf) const
    {
        const SharedSpace *SharedCLplusX = x->ObtainReadTempData("CLplusX");
        const double *CLplusX = SharedCLplusX->ObtainReadData();
        
        SharedSpace *SharediLplusX = new SharedSpace(3, n, n, num); // used to store the inverse of L
        double *iLplusX = SharediLplusX->ObtainWriteEntireData();
        
        
        const SharedSpace *SharedXplusCL = x->ObtainReadTempData("XplusCL");
        const double *XplusCL = SharedXplusCL->ObtainReadData();
        
        SharedSpace *SharediXplusL = new SharedSpace(3, n, n, num); // used to store the inverse of L
        double *iXplusL = SharediXplusL->ObtainWriteEntireData();

        
        const double *xM = x->ObtainReadData();
        double *gfVT = gf->ObtainWriteEntireData();
        for (integer i = 0; i < n * n; i++)gfVT[i] = 0.0;
        
        
        double *tmp = new double[n * n];
        double *LLplusX = new double[n * n * num];
        double *LXplusL = new double[n * n * num];
        
        double a1 = (1.0 - alpha)/2.0, a2 = (1.0 + alpha)/2.0;
        integer length = n * n, N = n, Nnum = n * n * num, info;
        double a, negone = -1.0;
        double *I = new double[n * n]; // make identity matrix
        for (integer i = 0; i < length; i++) I[i] = 0.0;
        for (integer i = 0; i < n; i++) I[i * n + i] = 1.0;

        
        const SharedSpace *Sharedcoef3 = x->ObtainReadTempData("coef3");
        const double *coef3 = Sharedcoef3->ObtainReadData();
        
        dcopy_(&Nnum, const_cast<double *>(CLplusX), &GLOBAL::IONE, LLplusX, &GLOBAL::IONE);
        dcopy_(&Nnum, const_cast<double *>(XplusCL), &GLOBAL::IONE, LXplusL, &GLOBAL::IONE);

        for (integer i = 0; i < num; i++)
        {
            // Compute inv(LplusX) -> iLplusX
            dcopy_(&length, I, &GLOBAL::IONE, iLplusX + n * n * i, &GLOBAL::IONE);
            dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, LLplusX + n * n * i, &N, iLplusX + n * n * i, &N, &info);
            
            // compute inv(LplusX) = inv(LplusX)^T * inv(LplusX)
            dcopy_(&length, iLplusX + n * n * i, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, tmp, &N, &GLOBAL::DZERO, iLplusX + n * n * i, &N);
            
            
            // Compute inv(XplusL) -> iXplusL
            dcopy_(&length, I, &GLOBAL::IONE, iXplusL + n * n * i, &GLOBAL::IONE);
            dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, LXplusL + n * n * i, &N, iXplusL + n * n * i, &N, &info);
            
            // compute inv(LplusX) = inv(LplusX)^T * inv(LplusX)
            dcopy_(&length, iXplusL + n * n * i, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, tmp, &N, &GLOBAL::DZERO, iXplusL + n * n * i, &N);
            
            // a2 * iLplusX + a1 * iXplusL -> tmp
            dcopy_(&length, iXplusL + n * n * i, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &a2, iLplusX + n * n * i, &N, I, &N, &a1, tmp, &N);
            
             a = 0.5/coef3[i];
             dscal_(&length, &a, tmp, &GLOBAL::IONE);
            
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, I, &N, &GLOBAL::DONE, gfVT, &N);
//             for(integer j = 0; j < n * n; j++) gfVT[j] = gfVT[j] + tmp[j];
            
        }
        
        // compute X * gfVT -> tmp
         dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *>(xM), &N, gfVT, &N, &GLOBAL::DZERO, tmp, &N);
            
        // compute tmp * X -> gfVT
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, const_cast<double *>(xM), &N, &GLOBAL::DZERO, gfVT, &N);

        // compute gfVT - num * X -> tmp
        dcopy_(&length, gfVT, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        dcopy_(&length, const_cast<double *>(xM), &GLOBAL::IONE, gfVT, &GLOBAL::IONE);
        
        double negnum = 0.0;
        for (integer i = 0; i < num; i++) negnum = negnum - 0.5/coef3[i];
        dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, I, &N, &negnum, gfVT, &N);

        
        x -> AddToTempData("iLplusX", SharediLplusX);
        x -> AddToTempData("iXplusL", SharediXplusL);
        delete[] I;
        delete[] tmp;
        delete[] LLplusX;
        delete[] LXplusL;
    };
    
    
    
    
    
//----------------------------------- Riemannian Hessian -----------------------------------
    void SPDMeanSymmLDOneParamL1::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const
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
    
    
//     void SPDMeanLDOneParamL1::RieHessianEig(Variable *x, double *EigVal) const
//     {
//         // compute the Hessian matrix H: H_ij = <e_i, Hessf(x)[e_j]>
//         // then compute the eigenvalue of H (symmetrize H finally!)
//         integer d = n * (n + 1)/2;
//         integer dim = d * d, N = n, length = n * n;
//         for (integer i = 0; i < d; i++) EigVal[i] = 0.0;
//
//         double *H = new double[d * d]; // Hessian matrix
//
//         SPDManifold *Mani = dynamic_cast<SPDManifold *> (Domain);
//         if (!x->TempDataExist("L"))
//         {
//             Mani->CholeskyRepresentation(x);
//         }
//         const SharedSpace *SharedL = x->ObtainReadTempData("L");
//         Variable *LElement = SharedL->GetSharedElement();
//         const double *L = LElement->ObtainReadData();
//
//
//         // define Basis
// 		double *e = new double[n * n];
//         double *E = new double[n * n * d];
//         double *B = new double[n * n * d]; // orthonormal basis under affine invariant metric
//         double *tmp = new double[n * n];
//
//         double *I = new double[n * n]; // make identity matrix
//         for (integer i = 0; i < length; i++) I[i] = 0.0;
//         for (integer i = 0; i < n; i++) I[i * n + i] = 1.0;
//
//
//         // produce the standard basis E first, then L * E * L' -> E
//         for (integer i = 0; i < n * n; i++) e[i] = 0;
//         for (integer i = 0; i < n; i++)
//         {
//             e[i + i * n] = 1;
// //             ForDebug::Print("e:", e + i * n, n, 1);
//         }
//
//         integer idx = 0;
//         for (integer i = 0; i < n; i++)
//         {
//             // ei * ei' -> Ei
//             dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &GLOBAL::IONE, &GLOBAL::DONE, e + i * n, &N, e + i * n, &N, &GLOBAL::DZERO, E + i * n * n, &N);
//
//             // L * Ei * L' -> Bi
//             // L * Ei -> Bi
//             dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *> (L), &N, E + i * n * n, &N, &GLOBAL::DZERO, B + i * n * n, &N);
//             // Bi * L' -> Bi
//             dcopy_(&length, B + i * n * n, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
//             dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, tmp, &N, const_cast<double *> (L), &N, &GLOBAL::DZERO, B + i * n * n, &N);
//
// //             ForDebug::Print("E:", E + i * n * n, n, n);
//             idx++;
//         }
//
//
// 		double r2 = 1.0/sqrt(2.0);
// 		for (integer i = 0; i < n; i++)
// 		{
// 			for (integer j = i + 1; j < n; j++)
// 			{
// 				// ei * ej' -> Eij
//                 dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &GLOBAL::IONE, &GLOBAL::DONE, e + i * n, &N, e + j * n, &N, &GLOBAL::DZERO, E + idx * n * n, &N);
//                 dcopy_(&length, E + idx * n * n, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
//
//                 // (tmp + tmp')/sqrt(2) -> Eij
//                 dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &r2, tmp, &N, I, &N, &r2, E + idx * n * n, &N);
//
//
//                 // L * Ei * L' -> Bi
//                 // L * Ei -> Bi
//                 dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *> (L), &N, E + idx * n * n, &N, &GLOBAL::DZERO, B + idx * n * n, &N);
//                 // Bi * L' -> Bi
//                 dcopy_(&length, B + idx * n * n, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
//                 dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, tmp, &N, const_cast<double *> (L), &N, &GLOBAL::DZERO, B + idx * n * n, &N);
//
// //              ForDebug::Print("E:", E + idx * n * n, n, n);
//                 idx++;
// 			}
// 		}
//
//
//
//         Vector *EMPTYEXTR = new SPDVector(n, n);
//         Vector *exBi = EMPTYEXTR->ConstructEmpty();
//         Vector *exBj = EMPTYEXTR->ConstructEmpty();
//         Vector *exHessBj = EMPTYEXTR->ConstructEmpty();
//
//
//         double *exBjM = exBj->ObtainWriteEntireData();
//         double *exBiM = exBi->ObtainWriteEntireData();
//
//
//         idx = 0;
//         for (integer i = 0; i < d; i++)
//         {
//             for (integer j = i; j < d; j++)
//             {
//                 // compute Hessf(x)[B_j]
//                 dcopy_(&length, B + i * n * n, &GLOBAL::IONE, exBiM, &GLOBAL::IONE); // <------
//                 dcopy_(&length, B + j * n * n, &GLOBAL::IONE, exBjM, &GLOBAL::IONE); // <------
//                 RieHessianEta(x, exBj, exHessBj);
//                 H[i + j * d] = Mani->Metric(x, exBi, exHessBj);
// //                 printf("idx: %d, %4f\n", idx, H[idx]);
//                 idx++;
//             }
//         }
//
//         for (integer i = 0; i < d; i++)
// 		{
// 			for (integer j = i + 1; j < d; j++)
// 			{
// 				H[j + i * d] = H[i + j * d];
// 			}
// 		}
//
//
//
//         ForDebug::Print("H1:", H, d, d);
//
//         // symmetrize the Hessian matrix
//         double onehalf = 0.5;
//         double *tmpd = new double[d * d];
//         double *Id = new double[d * d]; // make identity matrix
//         for (integer i = 0; i < dim; i++) Id[i] = 0.0;
//         for (integer i = 0; i < d; i++) Id[i * d + i] = 1.0;
//         dcopy_(&dim, H, &GLOBAL::IONE, tmpd, &GLOBAL::IONE);
//         dgemm_(GLOBAL::T, GLOBAL::N, &d, &d, &d, &onehalf, tmpd, &d, Id, &d, &onehalf, H, &d);
//
//
//         ForDebug::Print("Hsym:", H, d, d);
//
//         double *eigenvalues = new double[d + d * d];
//         double *eigenvectors = eigenvalues + d;
//         Matrix EE(eigenvalues, d, 1), VV(eigenvectors, d, d);
//         Matrix MMt(H, d, d);
//         Matrix::EigenSymmetricM(GLOBAL::L, MMt, EE, VV);
//
//         for (integer i = 0; i < d; i++)
//         {
//              printf("i: %d, I am eigenvalues: %.4f\n", i, eigenvalues[i]);
//             EigVal[i] = eigenvalues[i];
//         }
//
//         delete EMPTYEXTR;
//         delete[] e;
//         delete[] E;
//         delete[] I;
//         delete[] tmp;
//         delete[] tmpd;
//         delete[] eigenvalues;
//
//     };
    
    
}; /*end of ROPTLIB namespace*/
