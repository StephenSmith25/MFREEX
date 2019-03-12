#include "ShapeFunction/neighbours_materialpoint.h"

#include "matrix.h"
#include "matrix2.h"
#include <stdbool.h>


#ifndef TOLERRANCE 
#define TOLERRANCE 1e-12
#endif

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


IVEC * get_materialpoint_neighbours(IVEC * neighbours, MATERIAL_POINT * MP, MAT * nodes)
{

	int dim = nodes->n;
	int numnodes = nodes->m;
	double distance = 0;
	int count_neighbours = 0;
	double * x = MP->coords_n_1;


	if ( neighbours == IVNULL)
	{
		neighbours  = iv_get(10);
	}


	double xS[dim];

	for (int i = 0; i < numnodes; ++i)
	{


		if ( dim == 2)
		{
			xS[0] = x[0] - nodes->me[i][0]; 
			xS[1] = x[1] - nodes->me[i][1];

			distance = xS[0] * (MP->invMI->me[0][0] * xS[0] + MP->invMI->me[0][1]*xS[1]) + 
					   xS[1] * (MP->invMI->me[1][0] * xS[0] + MP->invMI->me[1][1]*xS[1]);
		}else{

			xS[0] = x[0] - nodes->me[i][0]; 
			xS[1] = x[1] - nodes->me[i][1];
			xS[2] = x[2] - nodes->me[i][2];

			distance = xS[0] * (MP->invMI->me[0][0] * xS[0] + MP->invMI->me[0][1]*xS[1] + MP->invMI->me[0][2]*xS[2]) + 
					   xS[1] * (MP->invMI->me[1][0] * xS[0] + MP->invMI->me[1][1]*xS[1] + MP->invMI->me[1][2]*xS[2]) + 
					   xS[2] * (MP->invMI->me[2][0] * xS[0] + MP->invMI->me[2][1]*xS[1] + MP->invMI->me[2][2]*xS[2]) ;
		}


		if ( distance <= 1.00)
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

		}

	return neighbours;

}

