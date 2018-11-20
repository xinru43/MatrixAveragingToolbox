
#include "Manifolds/SPDEucManifold/SPDEucVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDEucVector::SPDEucVector(integer row, integer col)
	{
		Element::Initialization(2, row, col);
	};

	SPDEucVector *SPDEucVector::ConstructEmpty(void) const
	{
		return new SPDEucVector(size[0], size[1]);
	};
}; /*end of ROPTLIB namespace*/
