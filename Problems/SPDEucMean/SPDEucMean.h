/*
This file defines the class for the problem: 
Computing the karcher mean of SPD manifold, i.e.,
min_{X \in SPD} sum dist^2(X, A_i)
where A_i \in SPD and dist(A, B) = \|logm(A^{-1/2} B A^{-1/2})\|_F.

Problem --> SPDMean

*/

#ifndef SPDEUCMEAN_H
#define SPDEUCMEAN_H

#include "Manifolds/SPDEucManifold/SPDEucManifold.h"
#include "Manifolds/SPDEucManifold/SPDEucVariable.h"
#include "Manifolds/SPDEucManifold/SPDEucVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Others/MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SPDEucMean : public Problem{
	public:
		/*The sample SPD matrices Ai are stored as their Cholesky matrix, i.e., Ai = Li Li^T.
		inLs is a n by n by num array of double numbers.*/
		SPDEucMean(double *inLs, integer inn, integer innum);
		virtual ~SPDEucMean();
		virtual double f(Variable *x) const;

		virtual void RieGrad(Variable *x, Vector *gf) const;

		/*Riemannian action of the Hessian has not been done yet*/
		virtual void RieHessianEta(Variable *x, Vector *etax, Vector *xix) const;

		/*The Euclidean gradient and Hessian are not done.*/
		//virtual void EucGrad(Variable *x, Vector *egf) const;
		//virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		double *Ls;
		integer n;
		integer num;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of STIEBROCKETT_H
