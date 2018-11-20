
#include "Problems/SPDEucMean/SPDEucMean.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDEucMean::SPDEucMean(double *inLs, integer inn, integer innum)
	{
		Ls = inLs;
		n = inn;
		num = innum;
	};

	SPDEucMean::~SPDEucMean(void)
	{
	};

	double SPDEucMean::f(Variable *x) const
	{
		SPDEucManifold *Mani = dynamic_cast<SPDEucManifold *> (Domain);
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
		for (integer i = 0; i < num; i++)
		{
			dcopy_(&length, const_cast<double *> (L), &GLOBAL::IONE, Ltmp, &GLOBAL::IONE);
			/*Solve the linear system Li X = Lx, i.e., X = Li^{-1} Lx. The solution X is stored in LiiLx.
			Note that Li is a lower triangular matrix.
			Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
			dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, Ls + n * n * i, &N, Ltmp, &N, &info);
			if (info != 0)
			{
				printf("The cholesky decompsotion in SPDEucMean::f failed with info:%d!\n", info);
			}
			dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, Ltmp, &N, Ltmp, &N, &GLOBAL::DZERO, logLXL + n * n * i, &N);
			Matrix MMt(logLXL + n * n * i, n, n);
			Matrix::LogSymmetricM(GLOBAL::L, MMt, MMt);
		}
		length = n * n * num;

		delete[] Ltmp;
		x->AddToTempData("logLXL", SharedlogLXL);

		double result = dnrm2_(&length, logLXL, &GLOBAL::IONE);
		result *= result;
		result /= 2.0 * num;
		return result;
	};

	void SPDEucMean::RieGrad(Variable *x, Vector *gf) const
	{
		const SharedSpace *SharedL = x->ObtainReadTempData("L");
		Variable *LElement = SharedL->GetSharedElement();
		const double *L = LElement->ObtainReadData();
        
        
        const SharedSpace *SharedlogLXL = x->ObtainReadTempData("logLXL");
		const double *logLXL = SharedlogLXL->ObtainReadData();
		double *gfVT = gf->ObtainWriteEntireData();
		for (integer i = 0; i < n * n; i++)
			gfVT[i] = 0;
		const double *xM = x->ObtainReadData();

		double *tmp = new double[n * n];
		integer N = n, length = n * n, info;
        
        double *I = new double[n * n]; // make identity matrix
        for (integer i = 0; i < length; i++) I[i] = 0.0;
        for (integer i = 0; i < n; i++) I[i * n + i] = 1.0;
        
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
				printf("The cholesky decompsotion in SPDEucMean::RieGrad failed with info:%d!\n", info);
			}
 			dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, I, &N, tmp, &N, &GLOBAL::DONE, gfVT, &N);
		}
        dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (L), &N, gfVT, &N, &info);
        dtrtrs_(GLOBAL::L, GLOBAL::T, GLOBAL::N, &N, &N, const_cast<double *> (L), &N, gfVT, &N, &info);
        
		double scalar = 1.0 / num;
		dscal_(&length, &scalar, gfVT, &GLOBAL::IONE);
		delete[] tmp;
        delete[] I;
	};

    
    
	void SPDEucMean::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const
	{
// 		printf("warning: SPDEucMean::RieHessianEta has not been implemented!\n");
		etax->CopyTo(xix);
	};
}; /*end of ROPTLIB namespace*/
