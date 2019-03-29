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



// n choose 3, n is atleast 3
static int nChoose3[13] = {
	1, // base case n = 3 
	4, // n = 4;
	10,
	20,
	35,
	56,
	84,
	120,
	165,
	220,
	286,
	364,
	455,
};



void find_combinations(int combinations[][3],IVEC * points, int num_points)
{
    int i, j, k;
    int count = 0;
    for (i = 0; i < num_points - 2; i++)
    {
        for (j = i + 1; j < num_points - 1; j++)
        {
            for (k = j + 1; k < num_points; k++)
            {
                combinations[count][0] =  points->ive[i];
                combinations[count][1] =  points->ive[j];
                combinations[count][2] =  points->ive[k];
			
				++count;
			}
        }

    }

    return ;
}

static inline double  crossProduct(double * p1 , double * p2)
{
	return p1[0]*p2[1] - p1[1]*p2[0];
}

static inline bool SameSide(double p1[2], double p2[2], double a[2], double b[2] )
{
	double b_a[2] = {b[0] - a[0], b[1] - a[1] };
	double p1_a[2] = {p1[0] - a[0], p1[1] - a[1]};
	double p2_a[2] = {p2[0] - a[0], p2[1] - a[1]};

	double res_1 = crossProduct(b_a, p1_a);
	double res_2 = crossProduct(b_a, p2_a);
	double dot_res = res_1*res_2;


	if ( dot_res >= 0 )
		return true;
	else
		return false;
}

static inline bool PointInTriangle (double * p, double * a, double * b, double * c)
{


	if (SameSide(p,a, b,c) && SameSide(p,b, a,c) && SameSide(p,c, a,b) )
        return true;
    else
		return false;
}


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
	

			MP->r_cutoff = MP->beta *  distances->ve[6];


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

	int min_num_neighbours = 6;

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


	MP->r_cutoff = MP->beta *distance[5];
	 MP->num_neighbours = neighbour_RangeSearch(MP->neighbours
	 		,cells, MP->coords_n_1, MP->r_cutoff, nodes);


	// // initial size of support domain 
	// double base_distance = distance[0];


	// MP->r_cutoff = MP->beta*base_distance;
	// MP->num_neighbours = neighbour_RangeSearch(MP->neighbours
	// 	,cells, MP->coords_n_1, MP->r_cutoff, nodes);

	// while ( MP->num_neighbours < min_num_neighbours)
	// {	
	// 	MP->r_cutoff = 1.2*MP->r_cutoff;
	// 	MP->num_neighbours = neighbour_RangeSearch(MP->neighbours
	// 	,cells, MP->coords_n_1, MP->r_cutoff, nodes);

	// }




	// // // Test conditions
	// bool volume_condition = false;
 // 	bool mat_point_contained = false;

 // 	// variables
 // 	double area_support = 0;

	// while ( volume_condition == false){

	// 	area_support = MP->r_cutoff * MP->r_cutoff * PI;

	// 	if ( area_support >= MP->volume){
	// 		volume_condition = true;
	// 	}else{
	// 	MP->r_cutoff = 1.1*MP->r_cutoff;
	// 	exit(0);


	// 	}

	// }

	// check particle disturbution is non-degenerate 



	// forms a non zero 

 // 	while ( mat_point_contained == false)
 // 	{


 // 		int index = MP->num_neighbours-3;

 // 		//printf("MP->num_neighbours = %d \n", MP->num_neighbours);
	//  	int num_combinations = nChoose3[MP->num_neighbours-3];

 //   		int combinations[num_combinations][3];

 //   		//printf("num_combinations = %d \n", num_combinations);

 // 		find_combinations(combinations,MP->neighbours, MP->num_neighbours);
	//  	i = 0;
	// 	while (( i < num_combinations) && (mat_point_contained == false))
	// 	{
	// 		int index_1 = combinations[i][0];
	// 		int index_2 = combinations[i][1];
	// 		int index_3 = combinations[i][2];

	// 		mat_point_contained = PointInTriangle(MP->coords_n_1, 
	// 			nodes->me[index_1], nodes->me[index_2],nodes->me[index_3]); 
	// 		++i;

	// 	}


	// 	if ( mat_point_contained == false)
	// 	{
		
	// 		//MP->beta = 1.05*MP->beta;
	// 		MP->r_cutoff = 1.2*MP->r_cutoff;
	// 		MP->num_neighbours = neighbour_RangeSearch(MP->neighbours
	// 		,cells, MP->coords_n_1, MP->r_cutoff, nodes);

	// 		//printf("num_neighbours = %d \n ", MP->num_neighbours);


	// 	}

	// }

	 // MP->r_cutoff = 1.2*MP->r_cutoff;
	 // MP->num_neighbours = neighbour_RangeSearch(MP->neighbours
	 // 		,cells, MP->coords_n_1, MP->r_cutoff, nodes);







	return 0;

}