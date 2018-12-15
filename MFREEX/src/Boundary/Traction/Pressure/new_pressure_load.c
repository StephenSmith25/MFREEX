#include "Boundary/Traction/Pressure/new_pressure_load.h"




pressure_boundary * new_pressure_load(IVEC * points,  meshfreeDomain * mFree)
{

	pressure_boundary * pB = malloc(1*sizeof(pressure_boundary));
	MAT * coords = mFree->nodes;
	int dim = coords->n;
	int num_points = points->max_dim;
	

	// points
	pB->points = points;
	pB->coords = m_get(num_points,dim);

	// set up traction coordinates
	for ( int i = 0 ; i < num_points ; i++)
	{
		pB->coords->me[i] = coords->me[pB->points->ive[i]];
	}


	// integration will be performed using a trapzoidal rule
	pB->sf_traction = mls_shapefunction(pB->coords, "linear", "cubic", 2, 1, mFree);


	// find segment weights and normals

	for ( int i = 0 ; i < num_points ; i++)
	{
		
	}



	return pB;
}


int free_pressure_boundary(pressure_boundary * pB)
{
	return 0;
}
