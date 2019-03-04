#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "matrix2.h"
#include <sys/time.h>
#include <omp.h>
#include "mat2csv.h"
#include "trigen.h"
#include "Integration/material_point.h"


#define QUADRATURE_ORDER 2

#define INTEGRATION_TYPE 1

int main(int argc, char const *argv[])
{


	// TRIANGULATE THE DOMAIN
	char opt[20] = "pDq0a0.05";
	char fileName[30] = "square";
	double * points_out ;
	int * boundaryNodes;
	int numBoundary;
	int * nodalMarkers;
	int  numnodes;

	struct triangulateio * tri  =  trigen(&points_out,&boundaryNodes,opt,fileName,&numnodes,&numBoundary,&nodalMarkers,NULL);

	double * tri_points = tri->pointlist;
	int * triangles = tri->trianglelist;
	int number_of_triangles = tri->numberoftriangles;
	int * quad_orders = malloc(number_of_triangles*sizeof(int));
	for ( int i = 0 ; i < number_of_triangles ; i++)
	{
		quad_orders[i] = QUADRATURE_ORDER;
	}	

	QUAD_TRIANGLE * quad_triangles = create_triangle_quadrature_points(tri_points, 
		triangles,number_of_triangles, quad_orders );


	printf("got here \n ");
	m_foutput(stdout, quad_triangles->QUAD_POINTS);



	int num_quad_points = quad_triangles->QUAD_POINTS->m;

	// export gauss points
	FILE * fp; 

	fp = fopen("quad_points.csv","w");
	for ( int i = 0 ; i < quad_triangles->QUAD_POINTS->m ; i++)
	{

		fprintf(fp,"%lf,%lf\n",quad_triangles->QUAD_POINTS->me[i][0],
			quad_triangles->QUAD_POINTS->me[i][1]);

	}

	fclose(fp);

	fp = fopen("nodes.csv","w");
	for ( int i =0 ; i < numnodes ; i++)
	{
		fprintf(fp, "%lf,%lf\n",tri_points[2*i],tri_points[2*i+1]);
	}
	fclose(fp);


	fp = fopen("triangles.csv","w");


	for ( int i = 0 ; i < number_of_triangles ; i++)
	{
		fprintf(fp,"%d,%d,%d\n",triangles[3*i],triangles[3*i+1],triangles[3*i+2]);
	}
	fclose(fp);



	printf("sum of volumes = %lf \n",v_sum(quad_triangles->VOLUMES));
	// export triangles
	MAT * xI = m_get(numnodes,2);


	// CREATE MATERIAL POINTS FROM SCNI OR TRIANGLE QUADRATURE


	/* code */
	return 0;
}



