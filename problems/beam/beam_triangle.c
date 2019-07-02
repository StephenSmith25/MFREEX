#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "matrix2.h"
#include "generate_voronoi.h"
#include "meshgrid.h"
#include <sys/time.h>
#include "Integration/SCNI/generate_scni.h"
#include "Integration/SCNI/generate_mscni.h"
#include "mls_shapefunction.h"
#include "setDomain.h"
#include "smoothstep.h"
#include "Force/Internal/internalForce_hyperelastic.h"
#include "Force/Internal/internalForce_hyperelastic_S.h"
#include "Boundary/Contact/contactDetection.h"
#include "Boundary/Displacement/setUpBC.h"
#include "Boundary/Displacement/enforceBC.h"
#include "Integration/SCNI/scni_update_B.h"
#include "Integration/material_point.h"
#include "Integration/mass_matrix.h"
#include "Integration/DomainMaterialPoint.h"
#include "mat2csv.h"
#include "trigen.h"



char * basis_type = "linear";
char * weight = "cubic";
char * kernel_shape = "radial";


int t_node = 1;

//char * integration_type = "TRIANGLE";
char * integration_type = "SCNI";




MATERIAL_TYPE material_type = HYPERELASTIC;

int is_AXI = 0;

// Meshfree parameters
	double dmax = 1.8;
const double dmax_x = 1.5;
const double dmax_y = 2.0;
double beta = 1.5;
int constant_support_size = 0;

double delta_t = 1e-6;


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
	double rho = 1000e-9;

	VEC * materialParameters = v_get(2);
	materialParameters->ve[0] = (E*nu)/((1+nu)*(1-2*nu));
	materialParameters->ve[1] =E/(2*(1+nu));


	// tip load
	double P = 10.00;
	// direction of the load
	VEC * dir_load = v_get(dim);
	dir_load->ve[0] = 0;
	dir_load->ve[1] = 1;

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


	// Read PLSG 
	char opt[20] = "pq20a0.1";
	char fileName[30] = "square";
	double * points_out ;
	int * boundaryNodes;
	int numBoundary;
	int * nodalMarkers;
	int  numnodes;
	// TRIANGULATION
	TRIANGLE * tri = trigen(opt,fileName);	

	boundaryNodes = tri->boundary;
	numBoundary = tri->num_boundary_points;
	nodalMarkers = tri->pointmarkers;
	numnodes = tri->num_points;
	points_out = tri->points;	


	MAT * xI = m_get(numnodes,dim);


	// create nodes
	for ( int i = 0 ; i < numnodes ; i++)
	{
		xI->me[i][0] = points_out[2*i];
		xI->me[i][1] = points_out[2*i+1];

	}

	struct timeval start, end;

	FILE * fp;
	fp = fopen("boundary.txt","w");
	for (int i = 0; i < numBoundary; ++i)
	{
		/* code */
		fprintf(fp,"%i\n",boundaryNodes[i]+1);
	}
	fclose(fp);


	/* ------------------------------------------*/
	/* ------------Meshfree Domain---------------*/
	/* ------------------------------------------*/

	// shape function parameters
	int compute = 3;
	VEC * dI = v_get(xI->m);

	double dmax_tensor[dim];
	dmax_tensor[0] = dmax_x;
	dmax_tensor[1] = dmax_y;


	// meshfree domain
	printf("got here \n");

	meshfreeDomain mfree = {.nodes = xI, .di = dI, .num_nodes = xI->m, .dim = dim, .IS_AXI = is_AXI,
		.weight_function = weight, .kernel_shape = kernel_shape, 
		.basis_type = basis_type,.is_constant_support_size = constant_support_size,
		.dmax_radial = dmax, .dmax_tensor = dmax_tensor, .beta = beta };
	
	setDomain(&mfree);



	// get transformation matrix at nodal points
	printf("getting shape functions \n ");

	shape_function_container * sf_nodes = mls_shapefunction(mfree.nodes,2, &mfree);

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






	// create saerch celsl
	BOUNDING_BOX * bounding_box = create_bounding_box(0, 21,
	0, 2, 0, 0);

	double cell_size[2] = {4,1.5};

	MAT * xI_copy = m_copy(xI,MNULL);
	CELLS * cells = create_cells(bounding_box, cell_size, dim, xI_copy);
	NODELIST * nodelist;

	FILE * fp1 = fopen("search_cells.csv","w");
	for ( int i = 0 ; i < cells->ny ; i++)
	{

		for ( int j = 0 ; j < cells->nx ; j++)
		{

		fprintf(fp1,"%lf,%lf,%lf,%lf\n",cells->cells[i][j].x[0],cells->cells[i][j].x[1],
			cells->cells[i][j].y[0],cells->cells[i][j].y[1]);
		
		}

	}
	fclose(fp1);

	fp1 = fopen("nodes.csv","w");
	for ( int i =0 ; i < numnodes ; i++)
	{
		for ( int k = 0 ; k < dim ; k++)
		{
			fprintf(fp1, "%lf",xI->me[i][k]);

			if ( k < dim - 1)
			{
				fprintf(fp1,",");
			}

		}
		fprintf(fp1, "\n");
	}
	fclose(fp1);




	/* ------------------------------------------*/
	/* ----------------MATERIAL POINTS----------------*/
	/* ------------------------------------------*/

	printf("-----------CREATING MATERIAL POINTS --------------\n");

	voronoi_diagram * vor = NULL;
	MATERIAL_POINTS * material_points = malloc(1*sizeof(material_points));


	int is_stabalised = 0;
	int is_AXI = 0;
	SCNI_OBJ * _scni_obj = NULL;

	int dim_s = dim;
	int dim_stress = dim*dim ;
	if ( is_AXI)
	{
		dim_s = dim_s +1;
		dim_stress = dim_stress + 1;
	}

	if ( strcmp ( integration_type, "SCNI") == 0)
	{
		material_points->num_material_points = numnodes;

		// generate clipped voronoi diagram
		vor = generate_voronoi(xI->base, boundaryNodes, xI->m, numBoundary, 2);
	 	// generate scni cells
		_scni_obj = generate_scni(vor, NULL , is_stabalised, is_AXI, dim, &mfree);

		// convert to material point representation
		int number_of_material_points = numnodes;

		material_points->MP = malloc(number_of_material_points*sizeof(MATERIAL_POINT*));
		material_points->num_material_points = number_of_material_points;


		for ( int i = 0 ; i < number_of_material_points ; i++)
		{
			material_points->MP[i] = malloc(1*sizeof(MATERIAL_POINT));
			material_points->MP[i]->inc_F = m_get(dim_s,dim_s);
			material_points->MP[i]->temp = m_get(dim_s,dim_s);

			material_points->MP[i]->temp_1 = m_get(dim_s,dim_s);
			material_points->MP[i]->index = i;
			material_points->MP[i]->fInt = _scni_obj->scni[i]->fInt;
			material_points->MP[i]->F_n = m_get(dim_s,dim_s);
			m_ident(material_points->MP[i]->inc_F);
			m_ident(material_points->MP[i]->F_n);
			material_points->MP[i]->Jn = 1;
			material_points->MP[i]->stressVoigt = v_get(dim_stress);
			material_points->MP[i]->coords_n = malloc(2*sizeof(double));
			material_points->MP[i]->coords_n_1 = malloc(2*sizeof(double));
			material_points->MP[i]->coords_n[0] = xI->me[i][0];
			material_points->MP[i]->coords_n[1] = xI->me[i][1];
			material_points->MP[i]->coords_n_1[0] = xI->me[i][0];
			material_points->MP[i]->coords_n_1[1] = xI->me[i][1];
			material_points->MP[i]->rho = rho;
			material_points->MP[i]->neighbours = _scni_obj->scni[i]->sfIndex;
			material_points->MP[i]->num_neighbours = _scni_obj->scni[i]->sfIndex->max_dim;
			material_points->MP[i]->shape_function = malloc(1*sizeof(shape_function));
			material_points->MP[i]->shape_function->phi = _scni_obj->scni[i]->phi;
			material_points->MP[i]->INTEGRATION_FACTOR = 1;
			material_points->MP[i]->volume = _scni_obj->scni[i]->area;
			material_points->MP[i]->stateNew = new_material_state(0, BUCKLEY,
					dim, is_AXI);
			material_points->MP[i]->stateOld = new_material_state(0, BUCKLEY,
					dim, is_AXI);
			material_points->MP[i]->B = _scni_obj->scni[i]->B;
		}

	
	}else if ( strcmp(integration_type, "TRIANGLE") == 0)
	{

	material_points = create_material_points(tri, 
			is_AXI, dim, integration_type, material_type,  cells, rho, beta,  &mfree);

	}


	int number_of_material_points = material_points->num_material_points;



	// form tempory Bmat

	// for ( int i = 0 ; i < numnodes ; i++)
	// {
		

	// 		// getting Bmat 
	// 		// Strain Displacement relationship 
	// 		material_points->MP[i]->B = BMAT(MNULL,
	// 			sf_nodes->sf_list[i],
	// 			dim,
	// 			0,
	// 			material_points->MP[i]->coords_n_1[0]);

	// 		material_points->MP[i]->neighbours = sf_nodes->sf_list[i]->neighbours;
	// 		material_points->MP[i]->num_neighbours =  sf_nodes->sf_list[i]->neighbours->max_dim;
	// 		material_points->MP[i]->fInt = v_resize(material_points->MP[i]->fInt, material_points->MP[i]->num_neighbours*2);


	// }

	/* ------------------------------------------*/
	/* -------------Boundary Conditions----------*/
	/* ------------------------------------------*/

	printf("Creating boundary conditions \n");
	// get boundary nodes
	IVEC * eb_nodes = iv_get(10);
	IVEC * traction_nodes = iv_get(10);

	int num_nodes_eb = 0;
	int num_nodes_trac = 0;
	double *x = NULL;
	for ( int i = 0 ; i < mfree.num_nodes ; i++ )
	{
		x = mfree.nodes->me[i];

		if ( x[0] < 1e-3)
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


	IVEC * temp = iv_copy(traction_nodes,IVNULL);
	iv_resize(traction_nodes, num_nodes_trac);
	traction_nodes->max_dim = num_nodes_trac;
	iv_move(temp, 0, num_nodes_trac, traction_nodes, 0);
	IV_FREE(temp);

	IVEC * temp_1 = iv_copy(eb_nodes,IVNULL);
	iv_resize(eb_nodes, num_nodes_eb);
	eb_nodes->max_dim = num_nodes_eb;
	iv_move(temp_1, 0, num_nodes_eb, eb_nodes, 0);
	IV_FREE(temp_1);


	iv_foutput(stdout,traction_nodes);

	iv_foutput(stdout, eb_nodes);





	// Traction nodes
	MAT * traction_nodes_coords = m_get(num_nodes_trac,dim);
	for ( int i = 0 ; i < num_nodes_trac ; i++)
	{
		traction_nodes_coords->me[i][0] = mfree.nodes->me[traction_nodes->ive[i]][0];
		traction_nodes_coords->me[i][1] = mfree.nodes->me[traction_nodes->ive[i]][1];

	}
	// get shape function and traction nodes 
	shape_function_container * phi_traction = mls_shapefunction(traction_nodes_coords,1, &mfree);

	/* ------------------------------------------*/
	/* ----------------Mass Vector---------------*/
	/* ------------------------------------------*/

	// VEC * nodal_mass = v_get(mfree.num_nodes);
	// VEC * inv_nodal_mass = v_get(mfree.num_nodes);

	// for ( int i = 0 ; i < mfree.num_nodes ; i++)
	// {
	// 	nodal_mass->ve[i] = _scni_obj->scni[i]->area*rho;
	// 	inv_nodal_mass->ve[i] = 1.000/nodal_mass->ve[i];

	// }

	VEC * nodal_mass = mass_vector(material_points, &mfree);

	VEC * inv_nodal_mass = v_get(mfree.num_nodes);

	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		inv_nodal_mass->ve[i] = 1.00/nodal_mass->ve[i];
	}


	printf("nodal mass = %lf \n", v_sum(nodal_mass));

	printf("----FINISHED BOUNDARY CONIDTION CREATION ---- \n");


	// /*  EB1  */
	EBC * eb1 = malloc(1*sizeof(EBC));
	eb1->dofFixed = 3;
	eb1->nodes = eb_nodes;
	int num_nodes_eb1 = eb1->nodes->max_dim;
	setUpBC(eb1,inv_nodal_mass,&mfree);


	/* ------------------------------------------*/
	/* ----------------Time stepping-------------*/
	/* ------------------------------------------*/
	
	// time parameters
	double t_max = 1; // 1s
	double t_n = 0;
	double t_n_1 = 0;

	// External Forces
	VEC * Fext_n_1 = v_get(dim*mfree.num_nodes);
	VEC * Fext_n = v_get(dim*mfree.num_nodes);

	// Internal Forces
	VEC * Fint_n_1 = v_get(dim*mfree.num_nodes);
	VEC * Fint_n = v_get(dim*mfree.num_nodes);

	// Net Force
	VEC * Fnet_n_1 = v_get(dim*mfree.num_nodes);
	eb1->uBar1 = v_get(eb1->nodes->max_dim);
	eb1->uBar2 = v_get(eb1->nodes->max_dim);

	// Kinematic variables
	// displacement
	VEC * d_n_1 = v_get(dim*mfree.num_nodes);
	VEC * d_n = v_get(dim*mfree.num_nodes);
	VEC * delta_disp = v_get(dim*mfree.num_nodes);
	VEC * nodal_disp = v_get(dim*mfree.num_nodes);
	// Velocity
	VEC * v_n_h = v_get(dim*mfree.num_nodes);
	VEC * v_n_mh = v_get(dim*mfree.num_nodes);
	VEC * v_n_1 = v_get(dim*mfree.num_nodes);
	VEC * v_n = v_get(dim*mfree.num_nodes);
	// Acceleration
	VEC * a_n_1 = v_get(dim*mfree.num_nodes);
	VEC * a_n = v_get(dim*mfree.num_nodes);

	MAT * updatedNodes = m_get(mfree.num_nodes,2);

	VEC * v_correct = v_get(dim*mfree.num_nodes);

	// tip load
	VEC * phi;
	IVEC * neighbours;
	double tipLoad = 0;

	// How often to write outputs
	int writeFreq = 1000;
	int fileCounter = 0;

	// time step counter;
	int n = 0;

	// Energy
	double Wext = 0;
	double Wint = 0;
	double Wkin = 0;
	double Wbal = 0;

	int num_dof = dim*mfree.num_nodes;



	MAT * nodes_X = m_copy(mfree.nodes,MNULL);



	struct timeval start3, end3;
	gettimeofday(&start3, NULL);
	double delta = 0;


	printf("-----STARTING TIME STEPPING ------\n");

	while ( t_n < t_max)
	//while ( n < 1)
	{

		// Update time step
		t_n_1 = t_n + delta_t;


		// Find d_n_1 and v_n_h
		// update velocity
		//v_mltadd(v_n_mh, a_n, delta_t, v_n_h);

		__mltadd__(v_n_h->ve, a_n->ve,delta_t,num_dof);

		// update displacements
		//v_mltadd(d_n,v_n_h,delta_t,d_n_1);
		__mltadd__(d_n_1->ve,v_n_h->ve,delta_t, num_dof);


		// // implement boundary conditions

		// for ( int i = 0 ; i < num_nodes_eb ; i++)
		// {
		// 	d_n_1->ve[eb_nodes->ive[i]*2] = 0;
		// 	d_n_1->ve[eb_nodes->ive[i]*2 + 1] = 0;
		// 	v_n_h->ve[eb_nodes->ive[i]*2] = 0;
		// 	v_n_h->ve[eb_nodes->ive[i]*2 + 1] = 0;
		// }

		v_zero(v_correct);

		/*  Implement BCs */
		enforceBC(eb1,d_n_1); 
		// find velocity correction
		sv_mlt(1.00/(delta_t),eb1->uCorrect1,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k] += v_correct->ve[k];
		}

		sv_mlt(1.000/(delta_t),eb1->uCorrect2,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k+1] += v_correct->ve[k];
		}


		// find new nodal positions
		mv_mlt(Lambda,d_n_1,nodal_disp);
		__add__(nodes_X->base, nodal_disp->ve, updatedNodes->base, num_dof);



		// Find the x and y points requried for plotting

		xPoint = nodal_disp->ve[traction_nodes->ive[t_node]*2+1]*xFactor;
		yPoint = tipLoad * pow(L,2)*(1/E)*(1/Ixx);


		/* ------------------------------------------*/
		/* ------------Find External Force-----------*/
		/* ------------------------------------------*/

		/*  Manually find point load */
		neighbours = phi_traction->sf_list[t_node]->neighbours;
		phi = phi_traction->sf_list[t_node]->phi;
		tipLoad = P*smoothstep(t_n_1,t_max,0);
		v_zero(Fext_n_1);
		for ( int i = 0 ; i < neighbours->max_dim; i++){
			// x
			Fext_n_1->ve[2*neighbours->ive[i]] = phi->ve[i]*tipLoad*dir_load->ve[0];
			 // y
			Fext_n_1->ve[2*neighbours->ive[i]+1] = phi->ve[i]*tipLoad*dir_load->ve[1];
		}

		/* ------------------------------------------*/
		/* ------------Find Internal Force-----------*/
		/* ------------------------------------------*/


		//double delta_t_min = internalForce_hyperelastic_S(Fint_n_1, _mscni_obj, d_n_1, v_n_h, materialParameters, "SVK", is_AXI, dim);
		// double delta_t_min = internalForce_hyperelastic(Fint_n_1, _scni_obj, d_n_1, v_n_h,
		// materialParameters, "SVK", is_AXI, dim,t_n_1);

		v_zero(Fint_n_1);
		for ( int i = 0 ; i < number_of_material_points ; i++)
		{

			internalForce_hyperelastic(Fint_n_1, material_points->MP[i], d_n_1,v_n_h,
 				materialParameters, &SVK , t_n_1);
		}



		/* ------------------------------------------*/
		/* ---------------Find Net Force-------------*/
		/* ------------------------------------------*/
		__sub__(Fext_n_1->ve, Fint_n_1->ve, Fnet_n_1->ve,num_dof);

		/* ------------------------------------------*/
		/* ---------------Find Acceleration----------*/
		/* ------------------------------------------*/

		// invM * Fne

		for ( int i = 0 ; i < numnodes  ; i++ )
		{
			a_n_1->ve[2*i] = Fnet_n_1->ve[2*i]*inv_nodal_mass->ve[i];
			a_n_1->ve[2*i+1] = Fnet_n_1->ve[2*i+1]*inv_nodal_mass->ve[i];
		}


		// update to interger time step velocities
		v_mltadd(v_n_h,a_n_1,0.5*delta_t,v_n_1);	


		for ( int i = 0 ; i < num_nodes_eb ; i++)
		{

			v_n_1->ve[eb_nodes->ive[i]*2] = 0;
			v_n_1->ve[eb_nodes->ive[i]*2 + 1] = 0;
		}


		// save outputs
		if ( n % writeFreq == 0 ){
			char filename[50];
			snprintf(filename, 50, "displacement_%d%s",fileCounter,".txt");
			mat2csv(updatedNodes,"./Displacement",filename);
			fp = fopen("loadDisp.txt","a");
			fprintf(fp,"%lf %lf\n",xPoint,yPoint);
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
		// v_copy(v_n_h,v_n_mh);
		v_copy(d_n_1,d_n);
		v_copy(v_n_1,v_n);
		v_copy(a_n_1,a_n);

		//delta_t = delta_t_min*0.85;

		t_n = t_n_1;

		++n;
		printf("%i  \t  %lf  \t %lf \t  %10.2E %10.2E \n",n,t_n,tipLoad, Wbal, delta_t);

	}

	gettimeofday(&end3, NULL);
	 delta = ((end3.tv_sec  - start3.tv_sec) * 1000000u + 
         end3.tv_usec - start3.tv_usec) / 1.e6;
	printf("Explicit routine took %lf seconds to run\n", delta);






	exit(0);


}
