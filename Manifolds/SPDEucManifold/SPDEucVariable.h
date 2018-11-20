/*
This file defines the class of a point on the the manifold of symmetric positive definite matrices (SPD)

SmartSpace --> Element --> SPDVariable

---- WH
*/

#ifndef SPDEUCVARIABLE_H
#define SPDEUCVARIABLE_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"
#include "Others/MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SPDEucVariable : public Element{
	public:
		/*Construct an empty variable on SPDEuc with only size information. */
		SPDEucVariable(integer n);

		/*Create an object of SPDEucVariable with same size as this SPDEucVariable.*/
		virtual SPDEucVariable *ConstructEmpty(void) const;

		/*This function randomly generates a point on the manifold.*/
		virtual void RandInManifold();
	};
}; /*end of ROPTLIB namespace*/
#endif // end of EUCVARIABLE_H
