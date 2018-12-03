#include "Integration/SCNI/scni_update_B.h"


static inline int isCCW(double * a, double * b, double * c){

	double det = (b[0]*c[1] + a[0]*b[1] + a[1]*c[0]) - (a[1]*b[0] + b[1]*c[0] + a[0]*c[1]);
	
	if ( det > 0 ){
		return 1;
	}else{
		return 0;
	}
}



// find if int a is in the list of integers stored in b
static inline int findInt( int a, int * b, int size_b)
{
	int i = 0;
	int b_temp = -1;

	for ( int i = 0 ; i < size_b ; i++)
	{
		int b_temp = b[i];


		if ( b_temp == a)
		{
			return i;
		}
	}

	return -1;
}

int scni_update_B(SCNI_OBJ * scni, VEC * disp, voronoi_diagram * voronoi, meshfreeDomain * Mfree){


	// cell points and meshfree nodes
	MAT * cell_verticies = voronoi->verticies;
	MAT * nodes = Mfree->nodes;

	int dim  = Mfree->dim;

	// get shape function at verticies
	shape_function_container * sf_verticies = mls_shapefunction(cell_verticies, "linear", "cubic", 2, 1, Mfree);





	return 0;
}
