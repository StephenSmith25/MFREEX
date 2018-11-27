#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "matrix2.h"
#include "generate_voronoi.h"
#include "meshgrid.h"
#include <sys/time.h>

// Jc_voronoi defintion
#define JC_VORONOI_IMPLEMENTATION
#define JCV_REAL_TYPE double
#define JCV_ATAN2 atan2
#define JCV_FLT_MAX 1.7976931348623157E+308
#include "jc_voronoi.h"
#include "Integration/SCNI/generate_scni.h"

#include "mls_shapefunction.h"
#include "setDomain.h"



int main(void )
{


	/* ------------------------------------------*/
	/* --------------Set up domain---------------*/
	/* ------------------------------------------*/
	
	// Dimension of problem
	int dim = 2;


	// create a 2D grid of test points on [0,1] x [0,1];
	int num_nodes_x = 16;
	int num_nodes_y = 16;

	// create nodes
	MAT * xI = meshgrid(0.000,10.000,num_nodes_x,0.000,10.000,num_nodes_y);
	int boundary[4];
	int count = 0;
	// find boundary nodes
	for ( int i = 0 ; i < xI->m; i++){
		double * x = xI->me[i];

		if ( (x[0] == 0) && ( x[1] == 0))
		{
			boundary[0] = i;
		}
		if ( (x[0] == 0) && ( x[1] == 10))
		{
			boundary[1] = i;
		}
		if ( (x[0] == 10) && ( x[1] == 10))
		{
			boundary[2] = i;
		}
		if ( (x[0] == 10) && ( x[1] == 0))
		{
			boundary[3] = i;
		}
	}
	struct timeval start, end;

	// generate clipped voronoi diagram
	gettimeofday(&start, NULL);
	voronoi_diagram * vor = NULL;
	vor = generate_voronoi(xI->base, boundary, xI->m, 4, 2);
 	// get time took to run
	gettimeofday(&end, NULL);
	double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6;

	// get time taken to run
	printf("voronoi took %lf seconds to run\n", delta);


	// Print voronoi diagram to file 
	FILE * fp = fopen("cells.txt","w");
	print_voronoi_diagram(fp,vor);
	fclose(fp);


	/* ------------------------------------------*/
	/* ------------Find Shape Functions----------*/
	/* ------------------------------------------*/

	// shape function parameters
	double dmax = 3;
	int constant_support_size = 1;
	char * basis = "quadratic";
	char * weight = "quartic";
	int compute = 3;
	VEC * dI = v_get(xI->m);

	// meshfree domain
	meshfreeDomain mfree = {.nodes = xI, .di = dI, .num_nodes = xI->m, .dim = dim};
	setDomain(&mfree,constant_support_size, dmax);

	// create a grid of sample points on [0,1] x [0,1]
	int num_sample_points_x = 50;
	int num_sample_points_y = 50;
	MAT * x_sample = meshgrid(0,10,num_sample_points_x,0,10,num_sample_points_y);


	struct timeval start1, end1;
	gettimeofday(&start1, NULL);
	// set up return container
	// generate shape functions at sample points 
	shape_function_container * sf_container = mls_shapefunction(x_sample, basis, weight, dim, compute, &mfree);
	gettimeofday(&end1, NULL);
	 delta = ((end1.tv_sec  - start1.tv_sec) * 1000000u + 
         end1.tv_usec - start1.tv_usec) / 1.e6;

	//free_shapefunction_container(sf_container);

	printf("Shape function took %lf seconds to run\n", delta);

	// v_foutput(stdout,sf_container->sf_list[0]->phi);
	// m_foutput(stdout,sf_container->sf_list[0]->dphi);
	// m_foutput(stdout,sf_container->sf_list[0]->d2phi);

	// V_FREE(dI);
	// M_FREE(x_sample);
	// M_FREE(xI);


	/* ------------------------------------------*/
	/* ----------------SCNI CELLS-----------------*/
	/* ------------------------------------------*/
	int is_stabalised = 0;
	int is_AXI = 0;
	SCNI ** _scni = NULL;


	gettimeofday(&start1, NULL);
	// set up return container
	// generate shape functions at sample points 
	_scni = generate_scni(vor, NULL , is_stabalised, is_AXI, dim, &mfree);
	gettimeofday(&end1, NULL);
	 delta = ((end1.tv_sec  - start1.tv_sec) * 1000000u + 
         end1.tv_usec - start1.tv_usec) / 1.e6;


	printf("scni function took %lf seconds to run\n", delta);


	V_FREE(dI);
	M_FREE(x_sample);
	M_FREE(xI);

	exit(0);


}
