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
#include "Integration/DomainMaterialPoint.h"
#include "cellSearch.h"
#include <stdbool.h>
#include <pthread.h>
#include "internal_force_mooney.h"
#include "input/read_input_mesh.h"
#include "input/read_config.h"
#include "Boundary/Displacement/DisplacementBC.h"



#ifndef DIM
#define DIM = 2
#endif 
/*  Function definitions */
int constant_support_size = 1;
char * basis_type = "linear";
char * weight = "cubic";
char * kernel_shape = "radial";


// how much larger can the domains get 
double beta =2;

// Meshfree parameters
const double dmax =4;
const double dmax_x = 1.5;
const double dmax_y = 1.5;
double tMax = 1;
double deltaT = 5e-8;

int writeFreq = 1000;
int printFreq = 1000;

const int dim = 2;
const int is_AXI = 0;

// Material type 
MATERIAL_TYPE material_type = HYPERELASTIC;
HYPERELASTIC_LAW hyperelastic_law = MOONEY_RIVLIN;
char * integration_type = "SCNI";



/*  Material  */

const double V_PRESCRIBED = 50;

const double T_RAMP_DISP =1e-2;

 
//#define IS_UPDATED
//#ifdef IS_UPDATED
//	#define UPDATE_FREQUENCEY 100
//#endif

double X_MIN_DIM = 0;
double X_MAX_DIM = 0.1;
double Y_MIN_DIM = 00.00;
double Y_MAX_DIM = 0.2;

int NUM_NODES_X = 15;
int NUM_NODES_Y = 15;



int main(int argc, char** argv) {



	/*////////////////////////////////////////////////////////// */
	/*                                                           */
	/*			        Initialisation                           */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */


	/* Initialsie random parameters  */
	FILE * fp;
	int i;

	int NUMBER_OF_THREADS = 3;
	if ( argv[2] != NULL)
	{
		NUMBER_OF_THREADS = atoi(argv[2]);

	}else{
		NUMBER_OF_THREADS = 1;	
	}


	/*////////////////////////////////////////////////////////// */
	/*                                                           */
	/*			        PREPRCOCESSING STAGE                     */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */

 
	/*  Material struct */
	VEC * materialParameters = v_get(3);
	materialParameters->ve[0] = 18.35;
	materialParameters->ve[1] = 1.468;
	materialParameters->ve[2] = 1.468e3;
	const double rho = 1.4089e-4;


	// CREATE POINTS BASED ON DELAUNAY TRIANGULATION
	char opt[20] = "p";
	char fileName[30] = "rubber";
	double * points_out ;
	int * boundaryNodes;
	int numBoundary;
	int * nodalMarkers;
	int  numnodes;
	double * temperatures;

	// TRIANGULATION
	TRIANGLE * tri = trigen(opt,fileName);	



	boundaryNodes = tri->boundary;
	numBoundary = tri->num_boundary_points;
	nodalMarkers = tri->pointmarkers;
	temperatures = tri->temperatures;
	numnodes = tri->num_points;
	points_out = tri->points;

	int * triangles = tri->triangles;
	int number_of_triangles = tri->num_triangles;


	fp = fopen("triangles.csv","w");
	for ( int i = 0 ; i < number_of_triangles ; i++)
	{
		fprintf(fp,"%d,%d,%d\n",triangles[3*i],triangles[3*i+1],triangles[3*i+2]);
		printf("%d,%d,%d\n",triangles[3*i],triangles[3*i+1],triangles[3*i+2]);
	}
	fclose(fp);

	MAT * xI = m_get(numnodes,dim);

	// create nodes
	for ( int i = 0 ; i < numnodes ; i++)
	{
		xI->me[i][0] = points_out[2*i];
		xI->me[i][1] = points_out[2*i+1];

	}

	printf("num_boundary = %d \n",numBoundary);
	fp = fopen("boundary.txt","w");
	for (int i = 0; i <numBoundary; ++i)
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



	// create stress points as barycenters of triangles

	MAT * stress_points = m_get(number_of_triangles,2);

	int number_of_stress_points = stress_points->m;


	//number_of_stress_points = 0;




	printf("number_of_stress_points %d \n",number_of_stress_points);


	for ( int i = 0 ; i < number_of_triangles ; i++)
	{
		int index_1 = triangles[3*i]-1;
		int index_2 = triangles[3*i+1]-1;
		int index_3 = triangles[3*i+2]-1;

		stress_points->me[i][0] = (xI->me[index_1][0] + xI->me[index_2][0] + xI->me[index_3][0])/3.00;
		stress_points->me[i][1] = (xI->me[index_1][1] + xI->me[index_2][1] + xI->me[index_3][1])/3.00;

		int indx = numnodes + i;

	}

	fp = fopen("stress_points.csv","w");
	for ( int i = 0 ; i < number_of_stress_points ; i++)
	{
		fprintf(fp,"%lf,%lf\n",stress_points->me[i][0],stress_points->me[i][1]);
	}
	fclose(fp);


	// all points together


	MAT * ALL_POINTS = m_get(numnodes + number_of_stress_points,2);
	for ( int i = 0 ; i < numnodes + number_of_stress_points ; i++)
	{

		if ( i < numnodes)
		{
			ALL_POINTS->me[i][0] = xI->me[i][0];
			ALL_POINTS->me[i][1] = xI->me[i][1];
		}else{
			int idx = i - numnodes;
			ALL_POINTS->me[i][0] = stress_points->me[idx][0];
			ALL_POINTS->me[i][1] = stress_points->me[idx][1];

		}

	}
	


	

	CELLS * cells = NULL;

	// CREATE NODES
	// 
	// BOUNDING_BOX * bounding_box = create_bounding_box(0, 3,
	// 0,3, 0, 0);

	// double cell_size[2] = {0.05,0.01};

	// MAT * xI_copy = m_copy(xI,MNULL);
	// CELLS * cells = create_cells(bounding_box, cell_size, dim, xI_copy);
	// NODELIST * nodelist;

	// fp = fopen("search_cells.csv","w");
	// for ( int i = 0 ; i < cells->ny ; i++)
	// {

	// 	for ( int j = 0 ; j < cells->nx ; j++)
	// 	{

	// 	fprintf(fp,"%lf,%lf,%lf,%lf\n",cells->cells[i][j].x[0],cells->cells[i][j].x[1],
	// 		cells->cells[i][j].y[0],cells->cells[i][j].y[1]);
		
	// 	}

	// }

	// fclose(fp);

	// int num_active_cells = 0;
	// active_cell * active_cells = get_active_cells(cells, &num_active_cells);
	// write_active_cells("active_cells.csv",active_cells);
	





	/* ------------------------------------------*/
	/* ------------Meshfree Domain---------------*/
	/* ------------------------------------------*/

	VEC * dI = v_get(xI->m);
	// meshfree domain
	meshfreeDomain mfree = {.nodes = xI, .di = dI, .num_nodes = xI->m, .dim = dim, .IS_AXI = is_AXI,
		.weight_function = weight, .kernel_shape = kernel_shape, 
		.basis_type = basis_type,.is_constant_support_size = constant_support_size,
		.dmax_radial = dmax};
	setDomain(&mfree);


	v_foutput(stdout,mfree.di);
	/* ---------------------------------------------------------------------------*/
	/* ---------------CREATE THE MATERIAL ( INTEGRATION) POINTS-------------------*/
	/* ---------------------------------------------------------------------------*/

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
		vor = generate_voronoi(ALL_POINTS->base, boundaryNodes,ALL_POINTS->m, numBoundary, 2);

		fp = fopen("cells.txt","w");
		print_voronoi_diagram(fp,vor);
		fclose(fp);



	 	// generate scni cells
		_scni_obj = generate_scni(vor, NULL , is_stabalised, is_AXI, dim, &mfree);

		// convert to material point representation
		int number_of_material_points = numnodes;

		material_points->MP = malloc(number_of_material_points*sizeof(MATERIAL_POINT*));
		material_points->num_material_points = number_of_material_points;
		material_points->dim = dim;
		material_points->IS_AXI = is_AXI;

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
			material_points->MP[i]->coords_n[0] = ALL_POINTS->me[i][0];
			material_points->MP[i]->coords_n[1] = ALL_POINTS->me[i][1];
			material_points->MP[i]->coords_n_1[0] = ALL_POINTS->me[i][0];
			material_points->MP[i]->coords_n_1[1] = ALL_POINTS->me[i][1];
			material_points->MP[i]->rho = rho;
			material_points->MP[i]->neighbours = _scni_obj->scni[i]->sfIndex;
			material_points->MP[i]->num_neighbours = _scni_obj->scni[i]->sfIndex->max_dim;
			material_points->MP[i]->shape_function = malloc(1*sizeof(shape_function));
			material_points->MP[i]->shape_function->phi = _scni_obj->scni[i]->phi;
			material_points->MP[i]->INTEGRATION_FACTOR = 1;
			material_points->MP[i]->volume = _scni_obj->scni[i]->area;
			material_points->MP[i]->stateNew = new_material_state(0, HYPERELASTIC,
					dim, is_AXI);
			material_points->MP[i]->stateOld = new_material_state(0, HYPERELASTIC,
					dim, is_AXI);
			material_points->MP[i]->B = _scni_obj->scni[i]->B;
		}











	
	}else if ( strcmp(integration_type, "TRIANGLE") == 0)
	{

		// material_points = create_material_points(tri, 
		// 		is_AXI, dim, integration_type, material_type,  cells, rho, beta,  &mfree);

	}


	int number_of_material_points = material_points->num_material_points;






	/* --------------------------------------------*/
	/* ----------------LUMPED MASSES---------------*/
	/* --------------------------------------------*/
	VEC * nodal_mass = mass_vector(material_points, &mfree);

	VEC * inv_nodal_mass = v_get(mfree.num_nodes);

	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		inv_nodal_mass->ve[i] = 1.00/nodal_mass->ve[i];
	}

	/* ------------------------------------------*/
	/* ----------------Boundaries---------------*/
	/* ------------------------------------------*/
	

	// Traction
	// int num_nodes_eb1=0;
	// int num_nodes_eb2=0;
	// int num_nodes_traction = 0;

	// for (int i = 0 ; i < mfree.num_nodes; i++)
	// {

	// 	double x = mfree.nodes->me[i][0];
	// 	double y = mfree.nodes->me[i][1];


	// 	if ( x < 1e-9)
	// 	{
	// 		++num_nodes_eb1;
	// 	}

	// 	if ( y < 1e-9)
	// 	{
	// 		++num_nodes_eb2;
	// 	}

	// 	double r = sqrt(x*x + y*y);
	// 	if ( r < 20+1e-9)
	// 	{
	// 		++num_nodes_traction;
	// 	}

	// }

	// IVEC * traction_nodes = iv_get(num_nodes_traction)

	/*  Set up essential boundary  */

	//m_foutput(stdout,pB->segment_normals);


	/*  Set up essential boundary  */
	EBC * eb1 = malloc(1*sizeof(EBC));
	EBC * eb2 = malloc(1*sizeof(EBC));


	int count1 = 0;
	int count = 0;

	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		double x = mfree.nodes->me[i][0];
		double y = mfree.nodes->me[i][1];


		if ( x < 1e-9)
		{
			++count;
		}
		if ( x > X_MAX_DIM - 1e-9)
		{
			++count1;
		}

	}

	printf("got here \n");

	printf("count = %d count1 = %d,\n",count,count1);

	eb1->nodes = iv_get(count);
	eb2->nodes = iv_get(count1);
	count1 = 0;
	count = 0;

	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		double x = mfree.nodes->me[i][0];
		double y = mfree.nodes->me[i][1];


		if ( x < 1e-9)
		{
			eb1->nodes->ive[count] = i;

			++count;
		}
		if ( x > X_MAX_DIM - 1e-9)
		{
			eb2->nodes->ive[count1] = i;
			++count1;
		}

	}
	eb1->dofFixed = 3;
	//getBoundary(&eb1->nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,3);
	//iv_addNode(eb1->nodes,traction_nodes->ive[traction_nodes->max_dim -1 ],'e');
	setUpBC(eb1,inv_nodal_mass,&mfree);

	m_foutput(stdout, eb1->coords);
	iv_foutput(stdout,eb1->nodes);


	/*  Set up essential boundary  */
	eb2->dofFixed = 3;
	//getBoundary(&eb2->nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,8);
	//iv_addNode(eb2->nodes,traction_nodes->ive[0],'s');
	setUpBC(eb2,inv_nodal_mass,&mfree);
	iv_foutput(stdout,eb2->nodes);
	m_foutput(stdout, eb2->coords);



	shape_function_container * sf_nodes = mls_shapefunction(mfree.nodes,1, &mfree);

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


	MAT * XI_n = m_get(numnodes,dim); 
	MAT * XI_n_1 = m_get(numnodes,dim);
	m_copy(mfree.nodes,XI_n);



	VEC * inc_disp = v_get(num_dof);

	VEC * D_N = v_get(num_dof);
	VEC * delta_disp = v_get(num_dof);
	VEC * nodal_disp = v_get(num_dof);

	/*  Force variables */
	VEC * Fext_n_1 = v_get(num_dof);
	VEC * Fext_n = v_get(num_dof);
	VEC * Fint_n_1 = v_get(num_dof);
	VEC * Fint_n = v_get(num_dof);
	VEC * Fnet = v_get(num_dof);
	VEC * v_correct = v_get(num_dof);
	eb1->uBar1 = v_get(eb1->nodes->max_dim);
	eb1->uBar2 = v_get(eb1->nodes->max_dim);
	eb2->uBar1 = v_get(eb2->nodes->max_dim);
	eb2->uBar2 = v_get(eb2->nodes->max_dim);




	VEC * R_pen = v_get(num_dof);

	MAT * nodes_X = m_copy(mfree.nodes,MNULL);

	// Energy
	double Wext = 0;
	double Wint = 0;
	double Wkin = 0;
	double Wbal = 0;
	double tStop = 0;


	double preStop = 0;


	double delta_t_min  = -1000;
	/*  Iteration counter */
	int n= 1;
	/*  File write counter */
	int fileCounter = 1;

	double pre_n_1 = 0; 


	VEC * FINT[NUMBER_OF_THREADS];
	VEC * RPEN[NUMBER_OF_THREADS];

	VEC * NODAL_MASS[NUMBER_OF_THREADS];
	for ( int i = 0 ; i < NUMBER_OF_THREADS ; i++)
	{

		NODAL_MASS[i] = v_get(numnodes);
		FINT[i] = v_get(num_dof);
		RPEN[i] = v_get(num_dof);

	}

	// damping
	double b1 = 0.06;
	double b2 = 1.44;
	double lambda = materialParameters->ve[0];
	double mu = materialParameters->ve[1];


	double u_pre = 0;
	double v_pre = 0;


	internal_force_args * INTERNAL_FORCE_ARGS = calloc(NUMBER_OF_THREADS,sizeof(internal_force_args));

	for ( int i = 0 ; i < NUMBER_OF_THREADS ; i++)
	{

		INTERNAL_FORCE_ARGS[i].NODAL_MASS = NODAL_MASS[i];
		INTERNAL_FORCE_ARGS[i].mass = nodal_mass;
		INTERNAL_FORCE_ARGS[i].FINT = FINT[i];
		INTERNAL_FORCE_ARGS[i].RPEN = RPEN[i];
		INTERNAL_FORCE_ARGS[i].material_points = material_points;
		INTERNAL_FORCE_ARGS[i].inc_disp = inc_disp;
		INTERNAL_FORCE_ARGS[i].materialParameters = materialParameters;
		INTERNAL_FORCE_ARGS[i].XI_n = XI_n;
		INTERNAL_FORCE_ARGS[i].XI_n_1 = XI_n_1;
		INTERNAL_FORCE_ARGS[i].cells = cells;
		INTERNAL_FORCE_ARGS[i].dt = deltaT;
		INTERNAL_FORCE_ARGS[i].velocity = v_n_h;


	}
	printf("NUMBER OF THREADS = %d \n", NUMBER_OF_THREADS);

	while ( t_n < tMax){


	/*  Update time step */
	t_n_1 = t_n + deltaT;

	// Find d_n_1 and v_n_h
	__mltadd__(v_n_h->ve, a_n_1->ve,deltaT,num_dof);
	// update displacements
	__mltadd__(d_n_1->ve,v_n_h->ve,deltaT, num_dof);
	



	// implement boundary conditions 
	// Apply boundary conditions as corrective accelerations
	v_pre = V_PRESCRIBED* smoothstep(t_n_1,T_RAMP_DISP,0);
	u_pre += v_pre*deltaT;

	if (u_pre > 10)
	{
		exit(0);
	}


	for ( int k = 0 ; k < eb2->nodes->max_dim ; k++)
	{
		eb2->uBar1->ve[k] = u_pre;
	}


	/*  Implement BCs */
	enforceBC(eb1,d_n_1); 
	// find velocity correction
	sv_mlt(1.00/(deltaT),eb1->uCorrect1,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k] += v_correct->ve[k];
		}

	sv_mlt(1.000/(deltaT),eb1->uCorrect2,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k+1] += v_correct->ve[k];
	}

	/*  Implement BCs */
	enforceBC(eb2,d_n_1); 
	//find velocity correction
	sv_mlt(1.00/(deltaT),eb2->uCorrect1,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k] += v_correct->ve[k];
	}

	sv_mlt(1.000/(deltaT),eb2->uCorrect2,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k+1] += v_correct->ve[k];
	}


	mv_mlt(Lambda,d_n_1,nodal_disp);
	__add__(nodes_X->base, nodal_disp->ve, XI_n_1->base, num_dof);



		//__add__(nodes_X->base, d_n_1->ve, XI_n_1->base, num_dof);

#ifdef IS_UPDATED
		move_nodes(cells, &active_cells, cell_size, XI_n_1);
#endif

	



		/* --------------------------------------------------------------*/
		/* --------------------------------------------------------------*/
		/*      			 INTERNAL FORCE ROUTINE 					 */
		/* --------------------------------------------------------------*/
		/* --------------------------------------------------------------*/
		__zero__(R_pen->ve,num_dof);
		__zero__(Fint_n_1->ve,Fint_n_1->max_dim);
		

		// update incremental displacemnt
		__sub__(d_n_1->ve, D_N->ve, inc_disp->ve, num_dof);


		for ( i = 0 ; i  < NUMBER_OF_THREADS ; i++)
		{
			__zero__(RPEN[i]->ve, num_dof);
			__zero__(FINT[i]->ve, num_dof);

		}
		#pragma omp parallel for num_threads(NUMBER_OF_THREADS) schedule(guided)
			for (i=0; i < number_of_material_points; i++){


				// assemble internal force
				int ID = omp_get_thread_num();
				INTERNAL_FORCE_ARGS[ID].MP = material_points->MP[i];
				internal_force_mooney(&INTERNAL_FORCE_ARGS[ID]);



				// update_material points
				#ifdef IS_UPDATED
				if ( n % UPDATE_FREQUENCEY == 0){
					material_points->MP[i] = update_material_point(material_points->MP[i], cells,
					XI_n_1, NODAL_MASS[ID]);
				}
				#endif


			}


		for (i=0; i < NUMBER_OF_THREADS; i++){
		
	 		__add__(RPEN[i]->ve, R_pen->ve,R_pen->ve, num_dof);
	 		__add__(FINT[i]->ve, Fint_n_1->ve,Fint_n_1->ve, num_dof);
		}


		if ( strcmp ( integration_type, "SCNI") == 0)
		{
			for ( int i = 0 ; i < numnodes ; i++)
			{
				material_points->MP[i]->coords_n_1[0] = XI_n_1->me[i][0];
				material_points->MP[i]->coords_n_1[1] = XI_n_1->me[i][1];	
			}
		}

		/*  Balance of forces */
		__sub__(Fint_n_1->ve, R_pen->ve, Fint_n_1->ve, num_dof);
		v_zero(Fext_n_1);
		__sub__(Fext_n_1->ve,Fint_n_1->ve,Fnet->ve,num_dof);

		// if ( n==10)
		// {
		// 	v_foutput(stdout, d_n_1);
		// 	printf("u_pre = %lf \n", u_pre);
		// 	exit(0);

		// }

		for (  i = 0 ; i < numnodes  ; i++ )
		{
			a_n_1->ve[2*i] = Fnet->ve[2*i]/nodal_mass->ve[i];
			a_n_1->ve[2*i+1] = Fnet->ve[2*i+1]/nodal_mass->ve[i];
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


		// used for updated rotuine
#ifdef IS_UPDATED
		if ( n % UPDATE_FREQUENCEY == 0){
			__zero__(nodal_mass->ve,numnodes);
			for ( i = 0 ; i  < NUMBER_OF_THREADS ; i++){
				__add__(NODAL_MASS[i]->ve, nodal_mass->ve,nodal_mass->ve, numnodes);
				__zero__(NODAL_MASS[i]->ve, numnodes);

			}
		}
#endif 

		// save outputs
		if (( n % writeFreq == 0 ) || (n == 0) ){

			char filename[50];
			snprintf(filename, 50, "displacement_%d%s",fileCounter,".txt");
			mat2csv(XI_n_1,"./Displacement",filename);
			

			snprintf(filename, 50, "MaterialPoints/materialpoints_%d%s",fileCounter,".txt");
			write_material_points(filename, material_points);

			// snprintf(filename, 50, "MaterialPoints/Domains/domains_%d%s",fileCounter,".txt");
			// write_domains(filename, material_points);	
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
		//deltaT = 0.85*delta_t_min; //delta_t_min;
		//v_copy(v_n_h,v_n_mh);
		v_copy(d_n_1,d_n);
		v_copy(v_n_1,v_n);

		#ifdef IS_UPDATED	
			if ( n % UPDATE_FREQUENCEY == 0){
				m_copy(XI_n_1,XI_n);
				v_copy(d_n_1,D_N);
		}
		#endif
		v_copy(a_n_1,a_n);
		t_n = t_n_1;
		// update iteration counter
		n++	;
		if ( n % printFreq == 0)
			printf("%i  \t  %lf  \t %lf \t  %10.2E %10.2E   %7lf     %7lf      %7lf \n",n,t_n,pre_n_1, Wbal,
			 deltaT,v_sum(Fext_n_1),v_sum(nodal_mass),v_sum(Fint_n_1));

		}

	/*////////////////////////////////////////////////////////// */
	/*////////////////////////////////////////////////////////// */





	/*  Program exit */
	exit(0);
}



	// #pragma omp parallel num_threads(NUMBER_OF_THREADS) 
	// 	{
	// 	// Internal forces
	// 	int ID = omp_get_thread_num();
	// 	IVEC * mat_points = material_point_groups[ID];
	// 	double P11,P12,P13,P21,P22,P23,P31,P32,P33;

	// 	MAT * F ;
	// 	MAT * Fdot;
	// 	VEC * stressVoigt ;
	// 	MAT * B ;
	// 	IVEC * neighbours;
	// 	MAT * F_r;
	// 	VEC * fInt;
	// 	double Xp_n, Xp_n_1, Yp_n, Yp_n_1;
	// 	int i = 0;

	// 	#pragma omp critical
	// 		__add__(NODAL_MASS[ID]->ve, nodal_mass->ve,nodal_mass->ve, numnodes);
		


	// 	__zero__(FINT[ID]->ve,num_dof);
	// 	__zero__(NODAL_MASS[ID]->ve,numnodes);
	// 	__zero__(RPEN[ID]->ve,num_dof);


		
	// 	//#pragma omp for nowait schedule(dynamic,4)
	// 	for ( int k = 0 ; k < mat_points->max_dim; k++)
	// 	{

	// 		i = mat_points->ive[k];

	// 		F = material_points->MP[i]->stateNew->F;
	// 		Fdot = material_points->MP[i]->stateNew->Fdot;
	// 		B = material_points->MP[i]->B;
	// 		neighbours = material_points->MP[i]->neighbours;
	// 		F_r = material_points->MP[i]->F_n;
	// 		stressVoigt = material_points->MP[i]->stressVoigt;
	// 		fInt = material_points->MP[i]->fInt;
	// 		int num_neighbours = material_points->MP[i]->num_neighbours;


	// 		// Compute incremental deformation gradient
	// 		defgrad_m(material_points->MP[i]->inc_F, B, neighbours, num_neighbours, inc_disp);
			

	// 		// compute total deformation gradient 
	// 		m_mlt(material_points->MP[i]->inc_F,
	// 			material_points->MP[i]->F_n,material_points->MP[i]->stateNew->F);
	
	
	// 		// Integration factor
	// 		double intFactor = material_points->MP[i]->volume*
	// 		material_points->MP[i]->INTEGRATION_FACTOR;


	// 		// Material law 
	// 		mooneyRivlin(stressVoigt,material_points->MP[i]->stateNew,materialParameters);
	// 		__zero__(fInt->ve,fInt->max_dim);


	// 		// push forward piola kirchoff stress to Omega_n configuration
	// 		sv_mlt(1.00/material_points->MP[i]->Jn,stressVoigt,stressVoigt);

	// 		P11 = stressVoigt->ve[0];
	// 		P22 = stressVoigt->ve[1];
	// 		P12 = stressVoigt->ve[2];
	// 		P21 = stressVoigt->ve[3];

	// 		stressVoigt->ve[0] = P11*F_r->me[0][0] + P12*F_r->me[0][1];
	// 		stressVoigt->ve[1] = P21*F_r->me[1][0] + P22*F_r->me[1][1];
	// 		stressVoigt->ve[2] = P11*F_r->me[1][0] + P12*F_r->me[1][1];
	// 		stressVoigt->ve[3] = P21*F_r->me[0][0] + P22*F_r->me[0][1];


	// 		vm_mlt(B,stressVoigt,fInt);



	// 		// Assemble internal force vector
	// 		for ( int k = 0 ; k < num_neighbours; k++){

	// 			int index = neighbours->ive[k];

	// 			FINT[ID]->ve[2*index] += intFactor * fInt->ve[2*k];
	// 			FINT[ID]->ve[2*index + 1] += intFactor *fInt->ve[2*k+1];


	// 		}


	// 		material_points->MP[i] = update_material_point(material_points->MP[i], cells,
	// 		 	XI_n_1, NODAL_MASS[ID]);

	// 		// Find error in displacment field
	// 		Xp_n_1 = material_points->MP[i]->coords_n_1[0];
	// 		Xp_n = material_points->MP[i]->coords_n[0];
	// 		Yp_n_1 = material_points->MP[i]->coords_n_1[1];
	// 		Yp_n = material_points->MP[i]->coords_n[1];


	// 		num_neighbours = material_points->MP[i]->num_neighbours;
	// 		neighbours = material_points->MP[i]->neighbours;




	// 		for ( int k = 0 ; k < num_neighbours; k++){

	// 			int index = neighbours->ive[k];

	// 			// Find vectors dX_n and dX_n_1
	// 			double delta_x_n[2] = {XI_n->me[index][0] - Xp_n,XI_n->me[index][1] - Yp_n};
	// 			double delta_x_n_1[2] = {XI_n_1->me[index][0] - Xp_n_1,XI_n_1->me[index][1] - Yp_n_1};

	// 			// Find error in displacement 
	// 			double x_tilde[2] = {material_points->MP[i]->inc_F->me[0][0]*(delta_x_n[0])
	// 				+ material_points->MP[i]->inc_F->me[0][1]*(delta_x_n[1]),
	// 				material_points->MP[i]->inc_F->me[1][0]*delta_x_n[0]
	// 				+ material_points->MP[i]->inc_F->me[1][1]*delta_x_n[1]};

	// 			double norm_x = sqrt(pow(delta_x_n[0],2) + pow(delta_x_n[1],2));
	// 			double e_x  = (delta_x_n_1[0] - x_tilde[0])/norm_x;
	// 			double e_y =  (delta_x_n_1[1] - x_tilde[1])/norm_x ;
		
	// 			// assemble forces 
	// 			RPEN[ID]->ve[2*index] += epsilon_penalty*material_points->MP[i]->shape_function->phi->ve[k]*e_x;
	// 			RPEN[ID]->ve[2*index+1] += epsilon_penalty*material_points->MP[i]->shape_function->phi->ve[k]*e_y;

	// 		}




	// 	} //end of loop over points
		

	// 	#pragma omp critical
	// 	{
	// 		__add__(RPEN[ID]->ve, R_pen->ve,R_pen->ve, num_dof);
	// 		__add__(FINT[ID]->ve, Fint_n_1->ve,Fint_n_1->ve, num_dof);
	// 	}


	// 	} // end of paralell region