#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "matrix2.h"
#include "generate_voronoi.h"
#include "meshgrid.h"
#include <sys/time.h>
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



	/* ------------------------------------------*/
	/* -------------Material & Loading-----------*/
	/* ------------------------------------------*/

	// material parameters and model (St. Venant Kirchoff)
	double nu = 0.0;
	double E = 4.8e3;
	char * material = "SVK";

	// tip load
	double P = 1.00;

	// Beam dimensions 
	double h = 1.00;
	double L = 20.00;

	// Plotting parameters
	double Ixx = 1.000/12.000;
	double yFactor = pow(L,2)/(E*Ixx);
	double xFactor = 1/L;
	double yPoint = P*yFactor;
	double xPoint = 0;


	/* ------------------------------------------*/
	/* --------------Tip Loaded Beam--------------*/
	/* ------------------------------------------*/

	// Number of nodes
	int num_nodes_x = 40;
	int num_nodes_y = 5;


	// create nodes
	printf("eneterd meshgrid \n ");
	MAT * xI = meshgrid(0.000,L,num_nodes_x,0.000,h,num_nodes_y);

	printf("got out of meshgrid \n ");
	int boundary[4];
	int count = 0;
	// find boundary nodes
	for ( int i = 0 ; i < xI->m; i++){
		double * x = xI->me[i];

		if ( (x[0] == 0) && ( x[1] == 0))
		{
			boundary[0] = i;
		}
		if ( (x[0] == 0) && ( x[1] == h))
		{
			boundary[1] = i;
		}
		if ( (x[0] == L) && ( x[1] == h))
		{
			boundary[2] = i;
		}
		if ( (x[0] == L) && ( x[1] == 0))
		{
			boundary[3] = i;
		}
	}

	printf("boundary = %d %d %d %d \n", boundary[0],boundary[1],boundary[2],boundary[3]);
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
	/* ------------Meshfree Domain---------------*/
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


	/* ------------------------------------------*/
	/* ----------------SCNI CELLS-----------------*/
	/* ------------------------------------------*/
	int is_stabalised = 0;
	int is_AXI = 0;
	SCNI ** _scni = NULL;
	struct timeval start2, end2;
	gettimeofday(&start2, NULL);
	// set up return container
	// generate shape functions at sample points 
	_scni = generate_scni(vor, NULL , is_stabalised, is_AXI, dim, &mfree);
	gettimeofday(&end2, NULL);
	 delta = ((end2.tv_sec  - start2.tv_sec) * 1000000u + 
         end2.tv_usec - start2.tv_usec) / 1.e6;
	printf("scni function took %lf seconds to run\n", delta);
	iv_foutput(stdout, _scni[2]->sfIndex);
	m_foutput(stdout,_scni[2]->B);



	// get boundary nodes
	IVEC * eb_nodes = iv_get(num_nodes_y);
	IVEC * traction_nodes = iv_get(num_nodes_y);

	int num_nodes_eb = 0;
	int num_nodes_trac = 0;
	double *x = NULL;
	for ( int i = 0 ; i < mfree.num_nodes ; i++ )
	{
		x = mfree.nodes->me[i];

		if ( x[0] == 0)
		{
			traction_nodes->ive[num_nodes_trac] = i;
			++num_nodes_trac;
		}
		if(  x[0] == L)
		{
			eb_nodes->ive[num_nodes_eb] = i;
			++num_nodes_eb;
		}

	}
	// 
	iv_foutput(stdout,traction_nodes);
	iv_foutput(stdout, eb_nodes);

	/* ------------------------------------------*/
	/* ----------------Time stepping-------------*/
	/* ------------------------------------------*/
	
	// time parameters
	double t_max = 1.00; // 1s
	double delta_t = 1e-6;
	double t_n = 0;
	double t_n_1 = 0;



	// External Forces
	VEC * Fext_n_1 = v_get(mfree.num_nodes);
	VEC * Fext_n = v_get(mfree.num_nodes);

	// Internal Forces
	VEC * Fint_n_1 = v_get(mfree.num_nodes);
	VEC * Fint_n = v_get(mfree.num_nodes);

	// Net Force
	VEC * Fnet_n_1 = v_get(mfree.num_nodes);


	// Kinematic variables
	// displacement
	VEC * d_n_1 = v_get(mfree.num_nodes);
	VEC * d_n = v_get(mfree.num_nodes);
	// Velocity
	VEC * v_n_h = v_get(mfree.num_nodes);
	VEC * v_n_mh = v_get(mfree.num_nodes);
	// Acceleration
	VEC * a_n_1 = v_get(mfree.num_nodes);
	VEC * a_n = v_get(mfree.num_nodes);


	while ( t_n < t_max)
	{

		// Update time step
		t_n_1 = t_n + delta_t;


		// Find d_n_1 and v_n_h

		// update velocity
		v_mltadd(v_n_mh, a_n, delta_t, v_n_h);
		// update displacements
		v_mltadd(d_n,v_n_h,delta_t,d_n_1);


		// Implement boundary conditions



		// Problem specific code
		// x-y 
		//xPoint = nodalDisp->ve[2*traction->nodes->ive[3]+1]*xFactor;
		//yPoint = traction->tMag * pow(L,2)*(1/Ey)*(1/Ixx);





	}









	exit(0);


}
