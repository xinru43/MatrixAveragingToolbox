/*
This file defines the class for the problem: 
Computing the karcher mean of SPD manifold, i.e.,
min_{X \in SPD} sum dist^2(X, A_i)
where A_i \in SPD and dist(A, B) = \|logm(A^{-1/2} B A^{-1/2})\|_F.

Problem --> SPDMean

---- WH
*/

#ifndef SPDMEANLDONEPARAM_H
#define SPDMEANLDONEPARAM_H

#include "SPDManifold.h"
#include "SPDVariable.h"
#include "SPDVector.h"
#include "Problem.h"
#include "SharedSpace.h"
#include "def.h"
#include "MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SPDMeanLDOneParam : public Problem{
	public:
		/*The sample SPD matrices Ai are stored as their Cholesky matrix, i.e., Ai = Li Li^T.
		inLs is a n by n by num array of double numbers.*/
		SPDMeanLDOneParam(double *inLs, integer inn, integer innum, double inalpha);
		virtual ~SPDMeanLDOneParam();
		virtual double f(Variable *x) const;

		virtual void RieGrad(Variable *x, Vector *gf) const;

		/*Riemannian action of the Hessian has not been done yet*/
		virtual void RieHessianEta(Variable *x, Vector *etax, Vector *xix) const;

		/*The Euclidean gradient and Hes sian are not done.*/
		//virtual void EucGrad(Variable *x, Vector *egf) const;
		//virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		double *Ls;
		integer n;
		integer num;
        double alpha;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of STIEBROCKETT_H
