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
#include "Integration/material_point.h"
#include "Integration/mass_matrix.h"

/*  Function definitions */

int constant_support_size = 0;
char * basis_type = "linear";
char * weight = "cubic";
char * kernel_shape = "radial";

// Meshfree parameters
const double dmax = 2.5;
const double dmax_x = 1.5;
const double dmax_y = 1.5;
double tMax = 1;

double deltaT = 5e-7;


const int dim = 2;
const int is_AXI = 0;
const double pressure = 230;




char * integration_type = "TRIANGLE";
char * MATERIAL_TYPE = "HYPERELASTIC";


/*  Material  */
char * material = "mooney_rivlin";
const double rho = 1000e-9;




int main(int argc, char** argv) {



	/*////////////////////////////////////////////////////////// */
	/*                                                           */
	/*			        Initialisation                           */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */


	/* Initialsie random parameters  */
	FILE * fp;
	int i;



	/*////////////////////////////////////////////////////////// */
	/*                                                           */
	/*			        PREPRCOCESSING STAGE                     */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */


	/*  Material struct */
	VEC * materialParameters = v_get(3);
	materialParameters->ve[0] = 1835.0;
	materialParameters->ve[1] = 146.8;
	materialParameters->ve[2] = 1e5;




	// CREATE POINTS BASED ON DELAUNAY TRIANGULATION


	char opt[20] = "pDYq0a0.1";
	char fileName[30] = "cylinder";
	double * points_out ;
	int * boundaryNodes;
	int numBoundary;
	int * nodalMarkers;
	int  numnodes;
	double * temperatures;


	// TRIANGULATION
	struct triangulateio * tri = trigen(&points_out,&boundaryNodes,opt,
		fileName,&numnodes,&numBoundary,&nodalMarkers,&temperatures);	


	MAT * xI = m_get(numnodes,dim);
	// create nodes
	for ( int i = 0 ; i < numnodes ; i++)
	{
		xI->me[i][0] = points_out[2*i];
		xI->me[i][1] = points_out[2*i+1];

	}

	fp = fopen("boundary.txt","w");
	for (int i = 0; i < numBoundary; ++i)
	{
		/* code */
		fprintf(fp,"%i\n",boundaryNodes[i]+1);
	}
	fclose(fp);


	fp = fopen("nodes.csv","w");
	for ( int i =0 ; i < numnodes ; i++)
	{
		for ( int k = 0 ; k < dim ; k++)
		{
			fprintf(fp, "%lf",xI->me[i][k]);

			if ( k < dim - 1)
			{
				fprintf(fp,",");
			}

		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	/* ------------------------------------------*/
	/* ------------Meshfree Domain---------------*/
	/* ------------------------------------------*/

	VEC * dI = v_get(xI->m);
	double dmax_tensor[dim];

	dmax_tensor[0] = dmax_x;
	dmax_tensor[1] = dmax_y;

	// meshfree domain
	meshfreeDomain mfree = {.nodes = xI, .di = dI, .num_nodes = xI->m, .dim = dim, .IS_AXI = is_AXI,
		.weight_function = weight, .kernel_shape = kernel_shape, 
		.basis_type = basis_type,.is_constant_support_size = constant_support_size,
		.dmax_radial = dmax, .dmax_tensor = dmax_tensor};
	
	setDomain(&mfree);


	/* ---------------------------------------------------------------------------*/
	/* ---------------CREATE THE MATERIAL ( INTEGRATION) POINTS-------------------*/
	/* ---------------------------------------------------------------------------*/

	voronoi_diagram * vor = NULL;
	MATERIAL_POINTS * material_points = NULL;

	if ( strcmp ( integration_type, "SCNI") == 0)
	{

		// generate clipped voronoi diagram
		vor = generate_voronoi(xI->base, boundaryNodes, xI->m, numBoundary, 2);
		

		material_points = create_material_points(vor, 
			is_AXI, dim, integration_type, MATERIAL_TYPE, rho, &mfree);
		

		// Write material points to file
		write_material_points("materialpoints.csv", material_points);

		// Print voronoi diagram to file 
		fp = fopen("cells.txt","w");
		print_voronoi_diagram(fp,vor);
		fclose(fp);
	
	}else if ( strcmp(integration_type, "TRIANGLE") == 0)
	{

		material_points = create_material_points(tri, 
			is_AXI, dim, integration_type, MATERIAL_TYPE, rho,  &mfree);
		

		// Write material points to file
		write_material_points("materialpoints.csv", material_points);

		double * tri_points = tri->pointlist;
		int * triangles = tri->trianglelist;
		int number_of_triangles = tri->numberoftriangles;


		fp = fopen("triangles.csv","w");
		for ( int i = 0 ; i < number_of_triangles ; i++)
		{
			fprintf(fp,"%d,%d,%d\n",triangles[3*i],triangles[3*i+1],triangles[3*i+2]);
		}
		fclose(fp);


		printf("got here \n");




	}

	// Test divergence free condition
	IVEC * index ;
	MAT * B; 
	MAT * check_B = m_get(dim*dim,2);
	printf("checking scni \n \n");
	int checkPoint = 33;


	printf("checking divergence free condition at point %d\n", checkPoint);
	printf("with coordinates %lf %lf\n", mfree.nodes->me[checkPoint][0], mfree.nodes->me[checkPoint][1]);
	printf("cell area = %lf \n", material_points->MP[checkPoint]->volume);
	printf("and neighbours \n");
	for ( int i = 0 ; i < material_points->num_material_points ; i++)
	{
		index = material_points->MP[i]->neighbours;
		B = material_points->MP[i]->B;

		for (int k = 0 ; k < index->max_dim ; k++){
			int indx = index->ive[k];

			if ( indx == checkPoint)
			{
				check_B->me[0][0] += B->me[0][2*k]*material_points->MP[i]->volume;
				check_B->me[1][1] += B->me[1][2*k+1]*material_points->MP[i]->volume;
				check_B->me[2][0] += B->me[2][2*k]*material_points->MP[i]->volume;
				check_B->me[3][1] += B->me[3][2*k+1]*material_points->MP[i]->volume;
			}
		}
	}

	m_foutput(stdout, check_B);
	/* --------------------------------------------*/
	/* ----------------LUMPED MASSES---------------*/
	/* --------------------------------------------*/
	VEC * nodal_mass = mass_vector(material_points, &mfree);

	VEC * inv_nodal_mass = v_get(mfree.num_nodes);

	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		inv_nodal_mass->ve[i] = 1.00/nodal_mass->ve[i];
	}

	printf("total mass = %lf (g) \n", v_sum(nodal_mass)*1000);

	double area = 0;
	for ( int i = 0 ; i < material_points->num_material_points ; i++)
	{
		area += material_points->MP[i]->volume;
	}

	printf("total area = %lf \n", area);


	/* ------------------------------------------*/
	/* ----------------Boundaries---------------*/
	/* ------------------------------------------*/
	

	// Traction
	IVEC * traction_nodes ;
	getBoundary(&traction_nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,2);
	pressure_boundary * pB = new_pressure_boundary(traction_nodes, &mfree);
	m_foutput(stdout,pB->coords);
	pB->is_axi = is_AXI;
	/*  Set up essential boundary  */

	//m_foutput(stdout,pB->segment_normals);


	/*  Set up essential boundary  */
	EBC * eb1 = malloc(1*sizeof(EBC));
	eb1->dofFixed = 1;
	getBoundary(&eb1->nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,3);
	iv_addNode(eb1->nodes,traction_nodes->ive[traction_nodes->max_dim -1 ],'e');
	setUpBC(eb1,inv_nodal_mass,&mfree);

	/*  Set up essential boundary  */
	EBC * eb2 = malloc(1*sizeof(EBC));
	eb2->dofFixed = 2;
	getBoundary(&eb2->nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,8);
	iv_addNode(eb2->nodes,traction_nodes->ive[0],'s');
	setUpBC(eb2,inv_nodal_mass,&mfree);


	/*  Get transformation matrix */
	/* Lambda  */

	shape_function_container * sf_nodes = mls_shapefunction(mfree.nodes, 1, &mfree);

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

	free_shapefunction_container(sf_nodes);



	///////////////////////////////////////////////////////////////
	
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
	/*  Kinematic variables */
	VEC * a_n_1 = v_get(num_dof);
	VEC * a_n = v_get(num_dof);
	VEC * d_n_1 = v_get(num_dof);
	VEC * d_n = v_get(num_dof);
	VEC * v_n_1 = v_get(num_dof);
	VEC * v_n_h = v_get(num_dof);

	VEC * v_n_mh = v_get(num_dof);
	VEC * v_n = v_get(num_dof);
	MAT * updatedNodes = m_get(num_dof,dim);


	VEC * delta_disp = v_get(num_dof);
	VEC * nodal_disp = v_get(num_dof);

	/*  Force variables */
	VEC * Fext_n_1 = v_get(num_dof);
	VEC * Fext_n = v_get(num_dof);
	VEC * Fint_n_1 = v_get(num_dof);
	VEC * Fint_n = v_get(num_dof);
	VEC * Fnet = v_get(num_dof);



	/* Boundary conditions */
	eb1->uBar1 = v_get(eb1->nodes->max_dim);
	eb1->uBar2 = v_get(eb1->nodes->max_dim);
	eb2->uBar1 = v_get(eb2->nodes->max_dim);
	eb2->uBar2 = v_get(eb2->nodes->max_dim);

	VEC * v_correct = v_get(num_dof);
	MAT * nodes_X = m_copy(mfree.nodes,MNULL);

	// Energy
	double Wext = 0;
	double Wint = 0;
	double Wkin = 0;
	double Wbal = 0;
	double tStop = 0;


	double preStop = 0;



	/*  Iteration counter */
	int n= 0;
	/*  File write counter */
	int writeFreq = 500;
	int fileCounter = 1;

	double pre_n_1 = 0; 

	
	/*  For writing to file */
	fp = fopen("loadDisp.txt","w");
	fprintf(fp,"%lf %lf\n",0.00,0.00);
	fclose(fp);

	while ( t_n < tMax){
	//while ( n < 10){
		/*  Update time step */
		t_n_1 = t_n + deltaT;

		// Find d_n_1 and v_n_h
		__mltadd__(v_n_h->ve, a_n_1->ve,deltaT,num_dof);
		// update displacements
		__mltadd__(d_n_1->ve,v_n_h->ve,deltaT, num_dof);
		// implement boundary conditions

			
		// Boundary conditions
		enforceBC(eb1,d_n_1); 
		// find velocity correction
		sv_mlt(1.00/(deltaT),eb1->uCorrect1,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k] += v_correct->ve[k];
		}

		enforceBC(eb2,d_n_1); 
		sv_mlt(1.000/(deltaT),eb2->uCorrect2,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k+1] += v_correct->ve[k];
		}


		// find new nodal positions
		mv_mlt(Lambda,d_n_1,nodal_disp);
		__add__(nodes_X->base, nodal_disp->ve, updatedNodes->base, num_dof);

		//__mltadd__(d_n_1->ve,v_n_h->ve,deltaT, num_dof);


		/*  Update tip load */
		pre_n_1 = pressure * smoothstep(t_n_1,tMax,0);
		update_pressure_boundary(pB, updatedNodes);
		v_zero(Fext_n_1);

		assemble_pressure_load(Fext_n_1, pre_n_1, pB);
	
		double delta_t_min = internalForce_hyperelastic(Fint_n_1, material_points, d_n_1, v_n_h,
		 materialParameters, material, is_AXI, dim,t_n_1);

		
		/*  Balance of forces */
		v_sub(Fext_n_1,Fint_n_1,Fnet);
		
		for ( int i = 0 ; i < numnodes  ; i++ )
		{
			a_n_1->ve[2*i] = Fnet->ve[2*i]*inv_nodal_mass->ve[i];
			a_n_1->ve[2*i+1] = Fnet->ve[2*i+1]*inv_nodal_mass->ve[i];
		}


		// update to interger time step velocities
		v_mltadd(v_n_h,a_n_1,0.5*deltaT,v_n_1);	


		for ( int i = 0 ; i < eb1->nodes->max_dim ; i++)
		{
			int indx = 2*eb1->nodes->ive[i];
			v_n_1->ve[indx] = 0;
		}


		for ( int i = 0 ; i < eb2->nodes->max_dim ; i++)
		{
			int indx = 2*eb2->nodes->ive[i]+1;
			v_n_1->ve[indx] = 0;
		}


		// save outputs
		if ( n % writeFreq == 0 ){
			char filename[50];
			snprintf(filename, 50, "displacement_%d%s",fileCounter,".txt");
			mat2csv(updatedNodes,"./Displacement",filename);

			fp = fopen("loadDisp.txt","a");
			printf("node = %d\n",traction_nodes->ive[traction_nodes->max_dim-1]) ;
			fprintf(fp,"%lf %lf\n",nodal_disp->ve[0],pre_n_1);		
			fclose(fp);
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
		deltaT = deltaT; //delta_t_min;
		v_copy(v_n_h,v_n_mh);
		v_copy(d_n_1,d_n);
		v_copy(v_n_1,v_n);
		v_copy(a_n_1,a_n);
		t_n = t_n_1;
		// update iteration counter
		n++	;

		printf("%i  \t  %lf  \t %lf \t  %10.2E %10.2E \n",n,t_n,pre_n_1, Wbal, deltaT);
	}


	/*////////////////////////////////////////////////////////// */
	/*////////////////////////////////////////////////////////// */





	/*  Program exit */
	return 0;
}
