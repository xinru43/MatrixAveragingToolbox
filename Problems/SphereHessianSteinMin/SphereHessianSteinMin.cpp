
#include "SphereHessianSteinMin.h"

/*Define the namespace*/
namespace ROPTLIB{

	SphereHessianSteinMin::SphereHessianSteinMin(double *inCs, double *inDs, integer inn, integer innum)
	{
        Cs = inCs;
        Ds = inDs;
        n = inn;
        num = innum;
	};

    SphereHessianSteinMin::~SphereHessianSteinMin(void)
	{
	};

double SphereHessianSteinMin::f(Variable *x) const
	{
        const double *xxM = x->ObtainReadData();
        SharedSpace *SharedCxD = new SharedSpace(3, n, n, num);
        double *CxD = SharedCxD->ObtainWriteEntireData();
        
        integer N = n, length = N * N;
        double *tmp = new double[n * n];
        double *tmp1 = new double[n * n];
        double *r = new double[num];
        
        for (integer i = 0; i < num; i++) r[i] = 0.0;
        
        
        for (integer i = 0; i < num; i++)
        {
//             printf("i: %d\n", i);
// //             ForDebug::Print("Cs:", Cs + n * n * i, n, n);
// //             ForDebug::Print("Ds:", Ds + n * n * i, n, n);
// //             ForDebug::Print("X:", const_cast<double *> (xxM), n, n);
//             
//             printf("Cs:\n");
//             for(integer j = 0; j < n * n; j++) printf("%.4f\n", Cs[j + n * n * i]);
//             
//             printf("Ds:\n");
//             for(integer j = 0; j < n * n; j++) printf("%.4f\n", Ds[j + n * n * i]);
//             
//             
//             printf("X:\n");
//             for(integer j = 0; j < n * n; j++) printf("%.4f\n", const_cast<double *> (xxM)[j]);

            
            // C * X -> CxD
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, Cs + n * n * i, &N, const_cast<double *> (xxM), &N, &GLOBAL::DZERO, CxD + n * n * i, &N);

            
            // CxD * D -> CxD
            dcopy_(&length, CxD + n * n * i, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, Ds + n * n * i, &N, &GLOBAL::DZERO, CxD + n * n * i, &N);
            
            // x*CxD -> tmp1
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *> (xxM), &N, CxD + n * n * i, &N, &GLOBAL::DZERO, tmp1, &N);

//             ForDebug::Print("xcxd:", tmp1, n, n);
            
            for (integer j = 0; j < n; j++) r[i] = r[i] + tmp1[j * n + j];
//             printf("trace: %.4f\n", r[i]);

            // output temp(:)^T * xxM(:), details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
//            r[i] = ddot_(&length, const_cast<double *> (xxM), &GLOBAL::IONE, CxD + n * n * i, &GLOBAL::IONE);
            
        }

        double result = 0;
        for (integer i = 0; i < num; i++) result = result + r[i];

        x->AddToTempData("CxD", SharedCxD);
        
        delete[] tmp;
        delete[] tmp1;
        delete[] r;
        
//         printf("f: %.4f\n", result);
		return result;
	};
    

	void SphereHessianSteinMin::EucGrad(Variable *x, Vector *egf) const
	{
        const double *xxM = x->ObtainReadData();

        const SharedSpace *SharedCxD = x->ObtainReadTempData("CxD");
        const double *CxD = SharedCxD->ObtainReadData();
        double *egfTV = egf->ObtainWriteEntireData();
        
        integer N = n, length = N * N;
        for (integer i = 0; i < length; i++) egfTV[i] = 0;
        
        double *I = new double[n * n]; // make identity matrix
        for (integer i = 0; i < length; i++) I[i] = 0.0;
        for (integer i = 0; i < n; i++) I[i * n + i] = 1.0;
        
        double scale = 2.0;
        
        double *DxC = new double[n * n];
        double *tmp = new double[n * n];
        
        
        for (integer i = 0; i < num; i++)
        {
//            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *> (CxD + n * n * i), &N, I, &N, &GLOBAL::DONE, egfTV, &N);
           
           
            // D * x -> DxC
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, Ds + n * n * i, &N, const_cast<double *> (xxM), &N, &GLOBAL::DZERO, DxC, &N);

            // DxC * C -> DxC
            dcopy_(&length, DxC, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, Cs + n * n * i, &N, &GLOBAL::DZERO, DxC, &N);
            
            // DxC + CxD -> tmp
            dcopy_(&length, const_cast<double *> (CxD + n * n * i), &GLOBAL::IONE, tmp, &GLOBAL::IONE);
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, DxC, &N, I, &N, &GLOBAL::DONE, tmp, &N);
            
            dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, tmp, &N, I, &N, &GLOBAL::DONE, egfTV, &N);
           
        }
        
        dcopy_(&length, egfTV, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
        double onehalf = 0.5;
        dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &onehalf, tmp, &N, I, &N, &onehalf, egfTV, &N);
        
        delete[] I;
        delete[] tmp;
        delete[] DxC;
        
	};

    
	void SphereHessianSteinMin::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
        std::cout << "warning: SphereHessianStein::EucHessianEta has not been implemented!" << std::endl;
        etax->CopyTo(exix);
	};
    


}; /*end of ROPTLIB namespace*/
