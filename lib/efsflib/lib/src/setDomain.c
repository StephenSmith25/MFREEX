#include "setDomain.h"




static inline double sq_distance(double *point_1, double * point_2, int  dim)
{
	double distance = 0;
	for (int i = 0; i < dim; ++i)
	{
		distance += pow(point_2[i] - point_1[i], 2);
		/* code */
	}

	return sqrt(distance);

}



int setDomain(meshfreeDomain * mfree, int constant_support_size, double dmax)
{

	MAT * nodes = mfree->nodes;


	int num_nodes = mfree->num_nodes;
	int dim = mfree->dim;



	double distance = 0;
	double distance_min = 1e6;
	// all nodes have same support size

	for (int i = 0; i < num_nodes; ++i)
	{
		distance_min = 1e6;
		distance = 1e6;


		// find distance from node I to each node J 
		for (int j = 0; j < num_nodes; ++j)
		{
			if ( i != j){
				// find distance to point;
				distance = sq_distance(nodes->me[i], nodes->me[j], dim);
				if (distance < distance_min )
				{
					distance_min = distance;
				}
			}
		}
		mfree->di->ve[i] = distance_min*dmax;
	}

	int index = 0;
	if  ( constant_support_size == 1)
	{

		double max_distance = v_max(mfree->di,&index);

		for (int i = 0; i < num_nodes; ++i)
		{
			mfree->di->ve[i] = max_distance;
			/* code */
		}

		return 0;



	// different support size for each node	
	}else{


		return 0;


	}

}


