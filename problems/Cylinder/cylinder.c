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
#include "thpool.h"



/*  Function definitions */

int constant_support_size = 0;
char * basis_type = "linear";
char * weight = "cubic";
char * kernel_shape = "radial";


// how much larger can the domains get 
const double alpha = 1.2;
double beta =1.2;

// Meshfree parameters
const double dmax =3;
const double dmax_x = 1.5;
const double dmax_y = 1.5;
double tMax = 0.4;
double deltaT = 5e-7;



int writeFreq = 1000;
int printFreq = 1000;

const int dim = 2;
const int is_AXI = 0;
const double pressure = 0.230;
const double T_RAMP_PRESSURE = 0.05;

// Material type 
MATERIAL_TYPE material_type = HYPERELASTIC;
HYPERELASTIC_LAW hyperelastic_law = MOONEY_RIVLIN;


char * integration_type = "TRIANGLE";


/*  Material  */
const double rho = 1000e-9;

 
#define IS_UPDATED
#ifdef IS_UPDATED
	#define UPDATE_FREQUENCEY 1000
#endif


int main(int argc, char** argv) {



	/*////////////////////////////////////////////////////////// */
	/*                                                           */
	/*			        Initialisation                           */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */


	/* Initialsie random parameters  */
	FILE * fp;
	int i;

	int NUMBER_OF_THREADS = 2;
	if ( argv[1] != NULL)
	{
		NUMBER_OF_THREADS = atoi(argv[1]);

	}else{
		NUMBER_OF_THREADS = 2;	
	}


	/*////////////////////////////////////////////////////////// */
	/*                                                           */
	/*			        PREPRCOCESSING STAGE                     */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */


	/*  Material struct */
	VEC * materialParameters = v_get(3);
	materialParameters->ve[0] = 1.8350;
	materialParameters->ve[1] = 0.1468;
	materialParameters->ve[2] = 1e3;






	// CREATE POINTS BASED ON DELAUNAY TRIANGULATION


	char opt[20] = "p";
	char fileName[30] = "cylinder";
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


	// CREATE NODES
	// 
	BOUNDING_BOX * bounding_box = create_bounding_box(0, 35,
	0, 35, 0, 0);

	double cell_size[2] = {1,1};

	MAT * xI_copy = m_copy(xI,MNULL);
	CELLS * cells = create_cells(bounding_box, cell_size, dim, xI_copy);
	NODELIST * nodelist;

	fp = fopen("search_cells.csv","w");
	for ( int i = 0 ; i < cells->ny ; i++)
	{

		for ( int j = 0 ; j < cells->nx ; j++)
		{

		fprintf(fp,"%lf,%lf,%lf,%lf\n",cells->cells[i][j].x[0],cells->cells[i][j].x[1],
			cells->cells[i][j].y[0],cells->cells[i][j].y[1]);
		
		}

	}

	fclose(fp);

	int num_active_cells = 0;
	active_cell * active_cells = get_active_cells(cells, &num_active_cells);
	write_active_cells("active_cells.csv",active_cells);
	
	IVEC * neighbours = iv_get(50);
	double x[2] = {13.5,16.12};
	int num_neighbours = neighbour_RangeSearch(neighbours, cells, x , 1, xI);
	iv_foutput(stdout, neighbours);
	printf("%d neighbours\n", num_neighbours);
	printf("Coordinates of neighbours \n");
	for ( int i = 0 ; i < num_neighbours ; i++)
	{
		int index = neighbours->ive[i];
		printf("%lf,%lf\n",xI->me[index][0], xI->me[index][1]);
	}







	// for ( int i = 0 ; i < numnodes; i++)
	// {
	// 	xI->me[i][0] += 11;
	// }

	//  printf("FIRST CALL : \n");
	// move_nodes(cells, &active_cells, cell_size, xI);
	
	// printf("SECOND CALL \n");
	// move_nodes(cells, &active_cells, cell_size, xI);

	// int count_2 = 1;
	// char filename_1[50];
	// char filename_2[50];




	// while ( count_2 < 50 ){


	// 	for ( int i = 0 ; i < numnodes; i++)
	// 	{
	// 		xI->me[i][0] += 0.2;
	// 		xI->me[i][1] += 0.2;

	// 	}

	// 	move_nodes(cells, &active_cells, cell_size, xI);


	// 	snprintf(filename_1, 50, "nodes_%d%s",count_2,".csv");
	// 	snprintf(filename_2, 50, "search_cells_%d%s",count_2,".csv");



	// 	fp = fopen(filename_1,"w");
	// 	for ( int i =0 ; i < numnodes ; i++)
	// 	{
	// 		for ( int k = 0 ; k < dim ; k++)
	// 		{
	// 			fprintf(fp, "%lf",xI->me[i][k]);

	// 			if ( k < dim - 1)
	// 			{
	// 				fprintf(fp,",");
	// 			}

	// 		}
	// 		fprintf(fp, "\n");
	// 	}
	// 	fclose(fp);

	// 	write_active_cells(filename_2,active_cells);

	// 	++count_2;


	// }



	// //snprintf(filename_2, 50, "nodes_%d%s",count_2-1,".csv");

	// // move_nodes(cells, &active_cells, cell_size, xI);
	// // move_nodes(cells, &active_cells, cell_size, xI);
	//  move_nodes(cells, &active_cells, cell_size, xI);
	// // move_nodes(cells, &active_cells, cell_size, xI);
	// // move_nodes(cells, &active_cells, cell_size, xI);
	// // move_nodes(cells, &active_cells, cell_size, xI);


	// write_active_cells(filename_2,active_cells);


	// fp = fopen(filename_1,"w");
	// 	for ( int i =0 ; i < numnodes ; i++)
	// 	{
	// 		for ( int k = 0 ; k < dim ; k++)
	// 		{
	// 			fprintf(fp, "%lf",xI->me[i][k]);

	// 			if ( k < dim - 1)
	// 			{
	// 				fprintf(fp,",");
	// 			}

	// 		}
	// 		fprintf(fp, "\n");
	// 	}
	// 	fclose(fp);

	// printf("num active cells = %d \n", num_active_cells);
	// // as points deform the active cell will change 
	// // so active cells should be a linked list to the cell index i,j
	// // and the next active cell along 

	// active_cell * activeCell = active_cells;

	// fp = fopen("active_nodes.csv","w");
	// NODELIST * current_p = NULL;
	// int nx = cells->nx;
	// int ny = cells->ny; 

	// while ( activeCell != NULL)
	// {	

	// 	int cell_number = activeCell->cell_number;
	// 	int i = (int) (floor(cell_number/nx));
	// 	int j = (int) (cell_number - i*(nx));
	// 	current_p = cells->cells[i][j].nodes;

	// 	while ( current_p != NULL)
	// 	{
	// 		fprintf(fp,"%d,",current_p->node_number+1);
	// 		current_p = current_p->next;



	// 	}

	// 	fprintf(fp,"\n");


	// 	activeCell = activeCell->next;
	// }
	// fclose(fp);





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
			is_AXI, dim, integration_type, material_type, cells, rho, beta,&mfree);
		

		// Write material points to file
		write_material_points("materialpoints.csv", material_points);

		// Print voronoi diagram to file 
		fp = fopen("cells.txt","w");
		print_voronoi_diagram(fp,vor);
		fclose(fp);
	
	}else if ( strcmp(integration_type, "TRIANGLE") == 0)
	{

		material_points = create_material_points(tri, 
			is_AXI, dim, integration_type, material_type,  cells, rho, beta,  &mfree);
		
		// Write material points to file
		write_material_points("materialpoints_0.csv", material_points);


		double * tri_points = tri->points;
		int * triangles = tri->triangles;
		int number_of_triangles = tri->num_triangles;


		fp = fopen("triangles.csv","w");
		for ( int i = 0 ; i < number_of_triangles ; i++)
		{
			fprintf(fp,"%d,%d,%d\n",triangles[3*i],triangles[3*i+1],triangles[3*i+2]);
		}
		fclose(fp);






	}



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
	IVEC * traction_nodes ;
	getBoundary(&traction_nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,2);
	pressure_boundary * pB = new_pressure_boundary(traction_nodes, &mfree);
	m_foutput(stdout,pB->coords);
	pB->is_axi = is_AXI;
	/*  Set up essential boundary  */

	//m_foutput(stdout,pB->segment_normals);


	/*  Set up essential boundary  */
	EBC * eb1 = malloc(1*sizeof(EBC));
	EBC * eb2 = malloc(1*sizeof(EBC));


	int count1 = 0;
	eb1->nodes = iv_get(4);
	eb2->nodes = iv_get(4);
	int count = 0;
	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		double x = mfree.nodes->me[i][0];
		double y = mfree.nodes->me[i][1];


		if ( x == 0)
		{
			eb1->nodes->ive[count] = i;

			++count;
		}
		if ( y == 0)
		{
			eb2->nodes->ive[count1] = i;
			++count1;
		}

	}
	eb1->dofFixed = 1;
	//getBoundary(&eb1->nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,3);
	//iv_addNode(eb1->nodes,traction_nodes->ive[traction_nodes->max_dim -1 ],'e');
	setUpBC(eb1,inv_nodal_mass,&mfree);

	m_foutput(stdout, eb1->coords);
	iv_foutput(stdout,eb1->nodes);


	/*  Set up essential boundary  */
	eb2->dofFixed = 2;
	//getBoundary(&eb2->nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,8);
	//iv_addNode(eb2->nodes,traction_nodes->ive[0],'s');
	setUpBC(eb2,inv_nodal_mass,&mfree);
	iv_foutput(stdout,eb2->nodes);
	m_foutput(stdout, eb2->coords);


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



	
	// v_zero(nodal_mass);

	// for (i = 0 ; i < material_points->num_material_points ; i++)
	// {

	// material_points->MP[i] = update_material_point(material_points->MP[i], cells, mfree.nodes, nodal_mass);


	// }


	// SPLIT MATERIAL POINTS
	int number_of_material_points = material_points->num_material_points;
	IVEC * material_point_groups[NUMBER_OF_THREADS];

	int number_per_group = number_of_material_points / NUMBER_OF_THREADS;
	int remainder = number_of_material_points % NUMBER_OF_THREADS;

	for ( int i = 0 ; i < NUMBER_OF_THREADS ; i++)
	{
		material_point_groups[i] = iv_get(number_per_group);

		if ( i == NUMBER_OF_THREADS -1 )
		{
			material_point_groups[i] = iv_get(number_per_group+remainder);
		}
	}

	count = 0;

	for ( int i = 0 ; i < NUMBER_OF_THREADS -1 ; i++  )
	{
		for ( int k = 0 ; k < number_per_group ; k++)
		{
			material_point_groups[i]->ive[k] = count;
			++count;
		}
	}
	for ( int k = 0 ; k < number_per_group + remainder ; k++)
	{
		material_point_groups[NUMBER_OF_THREADS-1]->ive[k] = count;
		++count;

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


	VEC * R_pen = v_get(num_dof);

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


	/*  For writing to file */
	fp = fopen("loadDisp.txt","w");
	fprintf(fp,"%lf %lf\n",0.00,0.00);
	fclose(fp);

	internal_force_args * INTERNAL_FORCE_ARGS = calloc(NUMBER_OF_THREADS,sizeof(internal_force_args));

	for ( int i = 0 ; i < NUMBER_OF_THREADS ; i++)
	{

		INTERNAL_FORCE_ARGS[i].NODAL_MASS = NODAL_MASS[i];
		INTERNAL_FORCE_ARGS[i].mass = nodal_mass;
		INTERNAL_FORCE_ARGS[i].FINT = FINT[i];
		INTERNAL_FORCE_ARGS[i].RPEN = RPEN[i];
		INTERNAL_FORCE_ARGS[i].mat_points = material_point_groups[i];
		INTERNAL_FORCE_ARGS[i].material_points = material_points;
		INTERNAL_FORCE_ARGS[i].inc_disp = inc_disp;
		INTERNAL_FORCE_ARGS[i].materialParameters = materialParameters;
		INTERNAL_FORCE_ARGS[i].XI_n = XI_n;
		INTERNAL_FORCE_ARGS[i].XI_n_1 = XI_n_1;
		INTERNAL_FORCE_ARGS[i].cells = cells;
		INTERNAL_FORCE_ARGS[i].dt = deltaT;
		INTERNAL_FORCE_ARGS[i].velocity = v_n_h;


	}

	// create thread pool
	//threadpool thpool = thpool_init(NUMBER_OF_THREADS);

 //    pthread_t threads[NUMBER_OF_THREADS];
 //   	pthread_attr_t attr[NUMBER_OF_THREADS];

 //    for ( int i = 0 ;i < NUMBER_OF_THREADS ; i++)
 //    {
 //       pthread_attr_init(&attr[i]);
 //       pthread_attr_setdetachstate(&attr[i], PTHREAD_CREATE_JOINABLE); 	
 //    }


	// threadpool thpool = thpool_init(NUMBER_OF_THREADS);
	while ( t_n < tMax){
	//for ( int v = 0 ; v < 20000 ; v++){


	/*  Update time step */
	t_n_1 = t_n + deltaT;

	// Find d_n_1 and v_n_h
	__mltadd__(v_n_h->ve, a_n_1->ve,deltaT,num_dof);
	// update displacements
	__mltadd__(d_n_1->ve,v_n_h->ve,deltaT, num_dof);
	// implement boundary conditions


	// Apply boundary conditions as corrective accelerations


		for ( int i =0 ; i < eb1->nodes->max_dim ; i++)
		{
			int index = eb1->nodes->ive[i];
			d_n_1->ve[2*index] = 0;
			v_n_h->ve[2*index] = 0;
		}

		for ( int i =0 ; i < eb2->nodes->max_dim ; i++)
		{
			int index = eb2->nodes->ive[i];
			d_n_1->ve[2*index+1] = 0;
			v_n_h->ve[2*index+1] = 0;
		}





		__add__(nodes_X->base, d_n_1->ve, XI_n_1->base, num_dof);

#ifdef IS_UPDATED
		move_nodes(cells, &active_cells, cell_size, XI_n_1);
#endif

		/* --------------------------------------------------------------*/
		/* --------------------------------------------------------------*/
		/*      				 EXTERNAL FORCE 					     */
		/* --------------------------------------------------------------*/
		/* --------------------------------------------------------------*/

		/*  Update pressureload */

	
		pre_n_1 = pressure* smoothstep(t_n_1,tMax,0);

		update_pressure_boundary(pB, XI_n_1);
		v_zero(Fext_n_1);
		assemble_pressure_load(Fext_n_1, pre_n_1, pB);



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

		/*  Balance of forces */
		__sub__(Fint_n_1->ve, R_pen->ve, Fint_n_1->ve, num_dof);
		__sub__(Fext_n_1->ve,Fint_n_1->ve,Fnet->ve,num_dof);


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
		if ( n % writeFreq == 0 ){

			char filename[50];
			snprintf(filename, 50, "displacement_%d%s",fileCounter,".txt");
			mat2csv(XI_n_1,"./Displacement",filename);
			write_material_points("materialpoints.csv", material_points);

			fp = fopen("loadDisp.txt","a");
			fprintf(fp,"%lf %lf\n",d_n_1->ve[0],pre_n_1);		
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
			printf("%i  \t  %lf  \t %lf \t  %10.2E %10.2E \n",n,t_n,pre_n_1, Wbal, deltaT);

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