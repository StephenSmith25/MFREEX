#include "stress_points.h"
static inline int isCCW(double * a, double * b, double * c){

	double det = (b[0]*c[1] + a[0]*b[1] + a[1]*c[0]) - (a[1]*b[0] + b[1]*c[0] + a[0]*c[1]);
	
	if ( det > 0 ){
		return 1;
	}else{
		return 0;
	}
}



STRESS_POINT ** create_stress_points(MAT * stress_point_coords,voronoi_diagram * voronoi, meshfreeDomain * Mfree)
{


	MAT * cell_verticies = voronoi->verticies;
	MAT * nodes = Mfree->nodes;



	int num_stress_points = stress_point_coords->m;
	STRESS_POINT ** stress_points = malloc(num_stress_points*sizeof(STRESS_POINT*));


	shape_function_container * sf_values = mls_shapefunction(stress_point_coords, 2, Mfree);

	double center[2];


	VEC * phi = VNULL;
	IVEC * neighbours = IVNULL;


	for ( int i = 0 ; i<num_stress_points ; i++)
	{
		int indx = Mfree->num_nodes + i;


		stress_points[i] = malloc(1*sizeof(STRESS_POINT));

		int * cell_index = voronoi->index[indx];
		int num_cell_verticies = voronoi->num_cell_verticies[indx];

		// re initialise parameters
		double area = 0;

		// check orientation of the polygon
		double * v1 = cell_verticies->me[cell_index[0]];
		double * v2 = cell_verticies->me[cell_index[1]];
		double * v3 = cell_verticies->me[cell_index[2]];

		// Get cell area 
		for ( int k = 0 ; k < num_cell_verticies ; k++)
		{
			v1 = cell_verticies->me[cell_index[k]];
			if ( k == (num_cell_verticies -1))
			{
				v2 = cell_verticies->me[cell_index[0]];

			}else{
				v2 = cell_verticies->me[cell_index[k+1]];

			}
			area += (v1[0]*v2[1] - v2[0]*v1[1]);
		}

		area = 0.5*fabs(area);


		// Generate Bmat
		stress_points[i]->B = BMAT(MNULL,sf_values->sf_list[i],2,1,stress_point_coords->me[i][0]);
		stress_points[i]->volume = area;
		stress_points[i]->r = stress_point_coords->me[i][0];
		stress_points[i]->fInt = v_get(2*sf_values->sf_list[i]->neighbours->max_dim);
		stress_points[i]->phi = v_copy(sf_values->sf_list[i]->phi,VNULL);
		stress_points[i]->neighbours = iv_copy(sf_values->sf_list[i]->neighbours,IVNULL);


	}
	free_shapefunction_container(sf_values);



	return stress_points;
}
