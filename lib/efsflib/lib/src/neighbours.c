#include "neighbours.h"

#include "matrix.h"
#include "matrix2.h"


double sq_distance(double *point_1, double * point_2, int  dim)
{
	double distance = 0;
	for (int i = 0; i < dim; ++i)
	{
		distance += pow(point_2[i] - point_1[i], 2);
		/* code */
	}

	return sqrt(distance);

}


IVEC * point_neighbours(double * x, meshfreeDomain * mfree)
{

	MAT * nodes = mfree->nodes;
	VEC * domainSize = mfree->di;

	
	int dim = mfree->dim;
	int numnodes = mfree->num_nodes;
	// approx number of neighbours is numnodes^(1/dim);
	IVEC * neighbours = iv_get((int)(floor ( sqrt(numnodes))));
	double distance = 0;
	int count_neighbours = 0;


	for (int i = 0; i < numnodes; ++i)
	{
		distance = sq_distance(x, nodes->me[i], dim);


		if ( distance < domainSize->ve[i])
		{
			neighbours->ive[count_neighbours] = i;
			count_neighbours = count_neighbours + 1;

			if ( count_neighbours >= neighbours->max_dim)
			{

				IVEC * temp = iv_get(neighbours->max_dim);
				iv_copy(neighbours,temp);
				iv_resize(neighbours,neighbours->max_dim+5);

				for (int k = 0 ; k < temp->max_dim; ++k)
				{
					neighbours->ive[k] = temp->ive[k];
				}
				IV_FREE(temp);

			}
		}
	}
	if ( neighbours->max_dim > count_neighbours)
	{
		IVEC * temp = iv_get(neighbours->max_dim);
		iv_copy(neighbours,temp);
		iv_resize(neighbours,count_neighbours);
		neighbours->max_dim = count_neighbours;

		for (int k = 0 ; k < count_neighbours; ++k)
		{
			neighbours->ive[k] = temp->ive[k];
		}
		IV_FREE(temp);
		return neighbours;

	}

	return neighbours;

}



