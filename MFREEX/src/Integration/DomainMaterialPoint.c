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
	

			MP->r_cutoff = 1.5* distances->ve[4];


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
	int i,j;


	int min_num_neighbours = 5; 

	double distance[MP->num_neighbours];

	// find min_num_neighbours for material point x_p
	for (j = 0; j < MP->num_neighbours; ++j)
	{
			// find distance to point;
			int index = MP->neighbours->ive[j];
			distance[j] = sq_distance(MP->coords_n_1, nodes->me[index], dim);

	}

	double a;
	for (i = 0 ; i < MP->num_neighbours ; i++)
	{
		for (  j = i+1 ; j < MP->num_neighbours ; j++)
		{
			if (distance[i] > distance[j])
			{
				a = distance[i];
				distance[i] = distance[j];
				distance[j] = a;
			}
		}
	}


	// printf("distance_min = %lf \n",distance_min);
	// if ( MP->num_neighbours < 6)
	// {
	// 	exit(0);
	// }
	MP->r_cutoff = 1.2*distance[4];
	//double area_support = MP->r_cutoff * MP->r_cutoff * PI;

	// if ( area_support < MP->volume)
	// {
	// 	printf("REACHED AREA CONDITION\n");
	// 	exit(0);
	// }

	MP->num_neighbours = neighbour_RangeSearch(MP->neighbours
		,cells, MP->coords_n_1, MP->r_cutoff, nodes);



	// while ( MP->num_neighbours < min_num_neighbours){


	// 	//MP->beta = 1.1*MP->beta;
	// 	MP->r_cutoff = 1.1*MP->r_cutoff;

	// 	MP->num_neighbours = neighbour_RangeSearch(MP->neighbours
	// 	,cells, MP->coords_n_1, MP->r_cutoff, nodes);

	// 	printf("got here \n");

	// }






	return 0;

}