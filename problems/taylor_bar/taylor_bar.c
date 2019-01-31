#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "matrix2.h"
#include "meshgrid.h"
#include "generate_voronoi.h"
#include "meshgrid.h"
#include <sys/time.h>
#include "Material/material.h"
#include "Integration/SCNI/generate_scni.h"
#include "Integration/SCNI/generate_mscni.h"
#include "mls_shapefunction.h"
#include "setDomain.h"
#include "smoothstep.h"
#include "Force/Internal/internalForce_hyperelastic.h"
#include "Boundary/Displacement/essential_boundary.h"
#include "mat2csv.h"
#include "trigen.h"
#include "Boundary/getBoundary.h"
#include "Boundary/iv_addNode.h"
#include "Boundary/Traction/Cavity/cavityVolume.h"
#include "Boundary/Traction/Cavity/flowRate.h"
#include "Boundary/Traction/Pressure/pressure_load.h"
#include <math.h>
#include "Deformation/poldec.h"
#include "Boundary/Contact/contactDetection.h"
#include "Boundary/Displacement/setUpBC.h"
#include "Boundary/Displacement/enforceBC.h"
#include "Integration/SCNI/scni_update_B.h"

// problem

const int dim = 2;
const int is_AXI = 1;

// Loading
const double INITIAL_BAR_VELOCITY = -300; // 300 m/s

// model
const double radius = 0.391/100; // 0.391cm
const double height = 2.346/100; // 2.346cm;
const double start_height = 0;
const double NUM_NODES_R = 10;
const double NUM_NODES_Z = 30;

// MATERIAL
const double YIELD_STRESS = 0.29e9; // Pa
const double mu = 0.3; // poisson ratio
const double E = 78.2e9; // Pa
const double rho = 2700; // kg/m^3

// time step
const double tMax = 10;
const double deltaT = 1e-9;

// Meshfree
const double dmax = 2;
const int constant_support_size = 1;
const int is_stabalised = 0;

// Define the rigid plate
const double RIGID_PLATE_WIDTH = 0.5/100;


// writing files
const int  writeFreq = 100;



int main(int argc, char** argv) {




	struct timeval beginPre, endPre, endExp;

	gettimeofday(&beginPre, NULL);



	/*////////////////////////////////////////////////////////// */
	/*////////////////////////////////////////////////////////// */




	char * material = "cubic_rivlin";

	VEC * materialParameters = v_get(4);
	materialParameters->ve[0] = 0.373e9;
	materialParameters->ve[1] = -0.031e9;
	materialParameters->ve[2] = 0.005e9;
	materialParameters->ve[3] = 1e9;


	/*////////////////////////////////////////////////////////// */
	/*                                                           */
	/*			        PREPRCOCESSING STAGE                     */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */

	// generate the model

	int num_boundary ;
	int * boundaryNodes;

	MAT * xI = meshgrid(0.00, radius, NUM_NODES_R, start_height, 
		height+start_height, NUM_NODES_Z, &boundaryNodes, &num_boundary);



	printf("creating voronoi diagram \n");
	// Construct voronoi diagram
	FILE * fp;
	struct timeval start, end;
	// generate clipped voronoi diagram
	gettimeofday(&start, NULL);
	voronoi_diagram * vor = NULL;
	vor = generate_voronoi(xI->base, boundaryNodes, xI->m, num_boundary, dim);
 	// get time took to run
	gettimeofday(&end, NULL);
	double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
		end.tv_usec - start.tv_usec) / 1.e6;

	// get time taken to run
	printf("voronoi took %lf seconds to run\n", delta);

	// Print voronoi diagram to file 
	fp = fopen("cells.txt","w");
	print_voronoi_diagram(fp,vor);
	fclose(fp);

	m_foutput(stdout, xI);
	// Rigid plate
	MAT * rigidPlate = m_get(2,2);
	rigidPlate->me[0][0] = 0;
	rigidPlate->me[1][0] = RIGID_PLATE_WIDTH;
	// rigidPlate->me[2][0] = 2*RIGID_PLATE_WIDTH/3.00;
	// rigidPlate->me[3][0] = RIGID_PLATE_WIDTH;


	/* ------------------------------------------*/
	/* ------------Meshfree Domain---------------*/
	/* ------------------------------------------*/

	// shape function parameters
	VEC * dI = v_get(xI->m);

	// meshfree domain
	meshfreeDomain mfree = {.nodes = xI, .di = dI, .num_nodes = xI->m, .dim = dim, .IS_AXI = is_AXI};
	setDomain(&mfree,constant_support_size, dmax);

	
	/* ------------------------------------------*/
	/* ------------------SCNI--------------------*/
	/* ------------------------------------------*/
	SCNI_OBJ * _scni_obj = NULL;
	struct timeval start2, end2;
	gettimeofday(&start2, NULL);
	_scni_obj = generate_scni(vor, NULL , is_stabalised, is_AXI, dim, &mfree);
	gettimeofday(&end2, NULL);
	delta = ((end2.tv_sec  - start2.tv_sec) * 1000000u + 
		end2.tv_usec - start2.tv_usec) / 1.e6;
	printf("scni function took %lf seconds to run\n", delta);


	/* ------------------------------------------*/
	/* ----------------Mass Vector---------------*/
	/* ------------------------------------------*/
	VEC * phi;
	IVEC * neighbours;
	VEC * nodal_mass = v_get(mfree.num_nodes);
	VEC * inv_nodal_mass = v_get(mfree.num_nodes);
	// get shape function and contact nodes
	shape_function_container * phi_nodes = mls_shapefunction(mfree.nodes, 
		"linear", "cubic", 2, 1, &mfree);

	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		phi = phi_nodes->sf_list[i]->phi;
		neighbours = phi_nodes->sf_list[i]->neighbours;

		double volume = _scni_obj->scni[i]->area;
		if ( is_AXI == 1)
		{
			double r = _scni_obj->scni[i]->r;
			volume = volume*r*2*PI;
		}
		for ( int k = 0 ; k < neighbours->max_dim ; k++)
		{
			int index = neighbours->ive[k];
			nodal_mass->ve[index] += phi->ve[k]*rho*volume;
		}

	}
	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		inv_nodal_mass->ve[i] = 1.00/nodal_mass->ve[i];
	}
	printf("total mass = %lf (g) \n", v_sum(nodal_mass)*1000);

	/* ------------------------------------------*/
	/* -----------Transformation Matrix----------*/
	/* ------------------------------------------*/

	shape_function_container * sf_nodes = mls_shapefunction(mfree.nodes, "linear", "cubic", 2, 1, &mfree);

	MAT * Lambda = m_get(2*mfree.num_nodes, 2*mfree.num_nodes);

	// u = Lambda * u_g
	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		VEC * phi = sf_nodes->sf_list[i]->phi;
		IVEC * neighbours  = sf_nodes->sf_list[i]->neighbours;
		for ( int k = 0 ; k < neighbours->max_dim ; k++)
		{
			Lambda->me[2*i][2*neighbours->ive[k]] += phi->ve[k]; 
			Lambda->me[2*i+1][2*neighbours->ive[k]+1] += phi->ve[k]; 

		}

	}


	// Find eb1 nodes
	IVEC * eb1_nodes = iv_get(NUM_NODES_Z);


	int count =0;
	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		double X = mfree.nodes->me[i][0];
		double Y = mfree.nodes->me[i][1];

		if ( X == 0)
		{
			eb1_nodes->ive[count] =  i;
			++count;
		}


	}

	iv_foutput(stdout, eb1_nodes);


	/*  Set up essential boundary  */
	EBC * eb1 = malloc(1*sizeof(EBC));
	eb1->nodes = eb1_nodes;
	eb1->dofFixed = 1;
	int num_nodes_eb1 = eb1->nodes->max_dim;
	setUpBC(eb1,inv_nodal_mass,&mfree);
	m_foutput(stdout,eb1->coords);

	/* ------------------------------------------*/
	/* ----------------Contact nodes------------*/
	/* ------------------------------------------*/


	/* ------------------------------------------*/
	/* ----------------Contact nodes------------*/
	/* ------------------------------------------*/

	IVEC * contact_nodes = iv_get(NUM_NODES_R);
	MAT * contact_nodes_coords = m_get(NUM_NODES_R,dim);
	count =0;
	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		double X = mfree.nodes->me[i][0];
		double Y = mfree.nodes->me[i][1];

		if ( Y == start_height)
		{
			contact_nodes->ive[count] =  i;
			contact_nodes_coords->me[count][0] = X;
			contact_nodes_coords->me[count][1] = Y;

			++count;
		}


	}


	/*  Set up essential boundary  */
	EBC * eb2 = malloc(1*sizeof(EBC));
	eb2->nodes = contact_nodes;
	eb2->dofFixed = 2;
	setUpBC(eb2,inv_nodal_mass,&mfree);



	printf("contact nodes = ");
	m_foutput(stdout, contact_nodes_coords);
		// get shape function and contact nodes
	shape_function_container * phi_contact = mls_shapefunction(contact_nodes_coords, 
		"linear", "quartic", 2, 1, &mfree);
	m_foutput(stdout,contact_nodes_coords);


	/*////////////////////////////////////////////////////////// */
	/*                                                           */
	/*			      EXPLICIT TIME STEPPNG                      */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */


	/*  time step  */
	double t_n_1;
	double t_n = 0;
	double t_n_h;

	int num_dof = mfree.num_nodes * dim;
	int numnodes = mfree.num_nodes;


	/*  Kinematic variables */
	VEC * a_n_1 = v_get(num_dof);
	VEC * a_n = v_get(num_dof);
	// Displacements
	VEC * d_n_1 = v_get(num_dof);
	VEC * d_n = v_get(num_dof);

	// Velocities
	VEC * v_n_1 = v_get(num_dof);
	VEC * v_n_h = v_get(num_dof);
	VEC * v_n_mh = v_get(num_dof);
	VEC * v_n = v_get(num_dof);


	// Displacmeent increment during the step
	VEC * delta_disp = v_get(num_dof);
	// Displacment of the nodes
	VEC * nodal_disp = v_get(num_dof);

	/*  Force variables */
	VEC * Fext_n_1 = v_get(num_dof);
	VEC * Fext_n = v_get(num_dof);
	VEC * Fint_n_1 = v_get(num_dof);
	VEC * Fint_n = v_get(num_dof);
	VEC * Fnet = v_get(num_dof);
	VEC * Fcont_n_1 = v_get(num_dof);
	VEC * Fcont_n = v_get(num_dof);

	for ( int i = 0 ; i < numnodes ; i++)
	{
		v_n_mh->ve[2*i+1 ] = INITIAL_BAR_VELOCITY;
	}



	VEC * v_correct = v_get(num_dof);


	// Original nodal postions
	MAT * nodes_X = m_copy(mfree.nodes,MNULL);
	// updated nodal positions
	MAT * updatedNodes = m_get(mfree.num_nodes,dim);

	// boundaries
	eb1->uBar1 = v_get(eb1->nodes->max_dim);
	eb2->uBar2 = v_get(eb2->nodes->max_dim);

	// Energy
	double Wext = 0;
	double Wint = 0;
	double Wkin = 0;
	double Wbal = 0;

	// Contact
	MAT * testPoint = m_get(1,dim);
	double distanceProj = 0.00;
	MAT * msNormal = m_get(1,dim);
	double f1Cor = 0.00;
	double f2Cor = 0.00;
	/*  Iteration counter */
	int n= 0;


	/*  File write counter */
		int fileCounter = 1;


	
	/*  For writing to file */
	fp = fopen("loadDisp.txt","w");
	fprintf(fp,"%lf %lf\n",0.00,0.00);
	fclose(fp);

	//while ( t_n < tMax){
	while ( n < 1000000){

		/*  Update time step */
		t_n_1 = t_n + deltaT;

		/*  Make a time step  */ 
		v_mltadd(v_n_mh,a_n,deltaT,v_n_h);
		v_mltadd(d_n,v_n_h,deltaT,d_n_1);


		__add__(nodes_X->base, d_n_1->ve, updatedNodes->base, num_dof);


		/*  Implement BCs */
		enforceBC(eb1,d_n_1); 
		// find velocity correction
		sv_mlt(1.00/(deltaT),eb1->uCorrect1,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k] += v_correct->ve[k];
		}

		// if ( n % 1000 == 0)
		// {
		// 	v_foutput(stdout, V_n)
		// }

		/*  Implement BCs */
		enforceBC(eb2,d_n_1); 
		// find velocity correction
		sv_mlt(1.00/(deltaT),eb2->uCorrect2,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			//v_n_h->ve[2*k+1] += v_correct->ve[k];
		}


		// find new nodal positions
		mv_mlt(Lambda,d_n_1,nodal_disp);
		__add__(nodes_X->base, nodal_disp->ve, updatedNodes->base, num_dof);


		// Find internal force

		double delta_t_min = internalForce_hyperelastic(Fint_n_1, _scni_obj, d_n_1, v_n_h,
		 materialParameters, "cubic_rivlin", is_AXI, dim,t_n_1);




		/*  Balance of forces */
		v_sub(Fext_n_1,Fint_n_1,Fnet);
		
		for ( int i = 0 ; i < numnodes  ; i++ )
		{
			a_n_1->ve[2*i] = Fnet->ve[2*i]*inv_nodal_mass->ve[i];
			a_n_1->ve[2*i+1] = Fnet->ve[2*i+1]*inv_nodal_mass->ve[i];
		}



		// update to interger time step velocities
		v_mltadd(v_n_h,a_n_1,0.5*deltaT,v_n_1);	


	


		// save outputs
		if ( n % writeFreq == 0 ){
			char filename[50];
			snprintf(filename, 50, "displacement_%d%s",fileCounter,".txt");
			mat2csv(updatedNodes,"./Displacement",filename);
			fileCounter++;	

		}
		/* ------------------------------------------*/
		/* -----------------Find Energy--------------*/
		/* ------------------------------------------*/


		// find change in displacement 
		__sub__(d_n_1->ve, d_n->ve, delta_disp->ve, num_dof);

		Wext += in_prod(delta_disp, Fext_n) + in_prod(delta_disp, Fext_n_1);
		Wint += in_prod(delta_disp, Fint_n) + in_prod(delta_disp, Fint_n_1);
		Wkin = 0;
		for ( int i = 0 ; i < mfree.num_nodes; i++ )
		{
			Wkin += 0.5*(v_n_1->ve[2*i]*nodal_mass->ve[i]*v_n_1->ve[2*i]);
			Wkin += 0.5*(v_n_1->ve[2*i+1]*nodal_mass->ve[i]*v_n_1->ve[2*i+1]);

		}

		double Wmax = max(Wext,Wint);
		Wmax =  max(Wmax,Wkin);
		if ( Wmax == 0)
		{
			Wbal = 0;
		}else{
			Wbal = fabs(Wkin + Wint - Wext)/Wmax;
		}
		// Store previous time step quanities for the kinematic, and force variables.
		v_copy(Fint_n_1,Fint_n);
		v_copy(Fext_n_1,Fext_n);
		v_copy(v_n_h,v_n_mh);
		v_copy(d_n_1,d_n);
		v_copy(v_n_1,v_n);
		v_copy(a_n_1,a_n);
		t_n = t_n_1;
		// update iteration counter
		n++	;

		printf("%i  \t  %lf  \t    %10.2E %10.2E \n",n,t_n, Wbal, deltaT);


	}


	/*////////////////////////////////////////////////////////// */
	/*////////////////////////////////////////////////////////// */



	/*////////////////////////////////////////////////////////// */
	/*			      POSTPROCESSING STAGE                       */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */

	
	/*  Free allocated memory */



	/*////////////////////////////////////////////////////////// */
	/*////////////////////////////////////////////////////////// */



	/*  Program exit */
	return 0;




}
