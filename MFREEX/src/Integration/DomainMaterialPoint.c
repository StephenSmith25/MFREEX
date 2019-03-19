/* ********************************
 * Author:       Stephen Smith 
 * License:	     Not yet decided 
 * Description:  Sets domain of influence for each material point
 *
 *//** @file DomainMaterialPoint.h *//*
 *
 ********************************/

#include "Integration/DomainMaterialPoint.h"
#include "ShapeFunction/neighbours_materialpoint.h"
#include "cellSearch.h"

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

int setDomainMaterialPoint(MAT * nodes, CELLS * cells, MATERIAL_POINT * MP)
{

	int num_nodes = nodes->m;
	int dim = nodes->n;



			VEC * distances = v_get(num_nodes); 
			PERM * order = px_get(num_nodes);
	


			// e.g find the 5 closest neighbours to x_p;
			for (int j = 0; j < num_nodes; ++j)
			{
				// find distance to point;
				distances->ve[j] = sq_distance(MP->coords_n_1, nodes->me[j], dim);


			}

			v_sort(distances,order);
	

			MP->r_cutoff = MP->beta * distances->ve[1];


			MP->num_neighbours = neighbour_RangeSearch(MP->neighbours,
			cells, MP->coords_n_1, MP->r_cutoff, nodes);



			PX_FREE(order);
			V_FREE(distances);

	


	return 0;

}

int updateDomainMaterialPoint(MAT * nodes, CELLS * cells,  MATERIAL_POINT * MP)
{

	int num_nodes = nodes->m;
	int dim = nodes->n;
	double distance_min = 1000;
	double distance = 0;


	// find min_num_neighbours for material point x_p
	for (int j = 0; j < MP->num_neighbours; ++j)
	{
			// find distance to point;
			int index = MP->neighbours->ive[j];



			distance = sq_distance(MP->coords_n_1, nodes->me[index], dim);

			if ( distance < distance_min)
				distance_min = distance;
	}


	MP->r_cutoff = MP->beta * distance_min;

	MP->num_neighbours = neighbour_RangeSearch(MP->neighbours
	,cells, MP->coords_n_1, MP->r_cutoff, nodes);





	return 0;

}