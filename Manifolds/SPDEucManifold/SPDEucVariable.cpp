
#include "Manifolds/SPDEucManifold/SPDEucVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDEucVariable::SPDEucVariable(integer n)
	{
		Element::Initialization(2, n, n);
	};

	SPDEucVariable *SPDEucVariable::ConstructEmpty(void) const
	{
		return new SPDEucVariable(size[0]);
	};

	void SPDEucVariable::RandInManifold(void)
	{
		integer n = size[0];
		/*temp is an n by n matrix*/
		double *temp = new double[n * n];
		for (integer i = 0; i < n; i++)
		{
			for (integer j = i; j < n; j++)
			{
				temp[i + j * n] = 0;
				temp[j + i * n] = genrandnormal();
			}
		}

		NewMemoryOnWrite();
		/*Space <-- temp * temp^T. Thus, Space points to a symmetric positive definite matrix. Therefore,
		this SPDEucVariable is a SPDEuc matrix.*/
		dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, temp, &n, temp, &n, &GLOBAL::DZERO, Space, &n);
		delete[] temp;
	};
}; /*end of ROPTLIB namespace*/
