#include "Boundary/Displacement/DisplacementBC.h"



void PrintConstraintType(DOF_CONSTRAINT * dof_constraint)
{
	DIRECTION dir = dof_constraint->dir;

	if ( dir == X)
	{
		printf("CONSTRAINED IN X DIRECTION\n");
	}
	if ( dir == Y)
	{
		printf("CONSTRAINED IN Y DIRECTION \n");
	}if ( dir == Z)
	{
		printf("CONSTRAINED IN Z DIRECTION \n");
	}



	return;
}


void ApplyBoundaryConditions();
