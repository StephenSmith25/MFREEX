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
#include "Force/Internal/internalForce_Inelastic_Buckley.h"
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
#include "PostProcess/saveDisp.h"
#include "input/read_input_points.h"
#include "Integration/mass_matrix.h"
#include "thpool.h"
#include "internal_force_buckley.h"



const int BUCKLEY_MATERIAL = 1;

// time step parameters
const double TMAX = 0.2;
double delta_t = 4e-7;

// Meshfree parameters
const double dmax = 2.5;
const double dmax_x =2;
const double dmax_y =2;
double beta = 1.4;


char * basis_type = "linear";
char * weight = "cubic";
char * kernel_shape = "radial";


const int is_stabalised = 0;
const int constant_support_size = 1;

// stretch rod
const double DISP_ROD_MAX = 90; // 132;


// Material type 
MATERIAL_TYPE material_type = BUCKLEY;
char * integration_type = "TRIANGLE";



//#define WITH_MOULD 
//define WITH_STRETCHROD


#define IS_UPDATED
#ifdef IS_UPDATED
	#define UPDATE_FREQUENCEY 1000
#endif

const int WRITE_FREQ =250;
const int PRINT_FREQ = 10;
int main(int argc, char** argv) {

	/*////////////////////////////////////////////////////////// */
	/*                                                           */
	/*			        Initialisation                           */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */



	int NUMBER_OF_THREADS = 2;
	if ( argv[1] != NULL)
	{
		NUMBER_OF_THREADS = atoi(argv[1]);

	}else{
		NUMBER_OF_THREADS = 2;	
	}



	FILE * fp;


	int dim = 2;
	int is_AXI = 1;

	/*  Material parameters */
	const double rho = 1380e-9;

	// material parameterss

	/*  Material struct */
	VEC * matParams = v_get(32);
	matParams->ve[0] = 2.814e-3; // VS
	matParams->ve[1] = 0.526e-3; // VP
	matParams->ve[2] = (1.7057e6); // mu*_0
	matParams->ve[3] = (328.76); // Tinf
	matParams->ve[4] = 358.15; // T*
	matParams->ve[5] = matParams->ve[4]; // Tf*
	matParams->ve[6] = (67.47); // Cv
	matParams->ve[7] = 1.23e5; // H0
	matParams->ve[8] = 8.314; // R
	matParams->ve[9] = 1.0e9; // Kb
	matParams->ve[10] = 6e8;// Gb
	// conformational constants
	matParams->ve[13] = 0.1553;// alpha_c
	matParams->ve[14] = 0.001;// eta_c
	matParams->ve[15] = 1.8098e17;// Ns_c
	matParams->ve[16] = 1.38e-17;// boltzmann constant kB
	// slippage
	matParams->ve[17] = 100;// lambdaCrit
	matParams->ve[18] = 383.15;// Ts 
	matParams->ve[19] = 0.653e6;// gamma0_ref = 0.653
	matParams->ve[20] = 10612;// Cs 10612
	matParams->ve[21] = 95.48;// Tinf 95.48
	matParams->ve[22] = 0.1565;// C1
	matParams->ve[23] = 39.937;// C2
	matParams->ve[24] = 0.9878;// beta
	matParams->ve[25] = 0.33;// poissons ratio

	// matParams->ve[17] = 100;// lambdaCrit
	// matParams->ve[18] = 383.15;// Ts 
	// matParams->ve[19] = 0.359e6;// gamma0_ref = 0.653
	// matParams->ve[20] = 7307.8;// Cs 10612
	// matParams->ve[21] = 152.95;// Tinf 95.48
	// matParams->ve[22] = 0.1565;// C1
	// matParams->ve[23] = 39.937;// C2
	// matParams->ve[24] = 0.9878;// beta
	// matParams->ve[25] = 0.33;// poissons ratio


	
	// crit lambda properties
	matParams->ve[26] = -0.0111; // C1
	matParams->ve[27] = 3.627; // C2
	matParams->ve[28] = 0.9856; // BETA
	matParams->ve[29] = -0.0356; // k
	matParams->ve[30] = 15.393; // b 
	matParams->ve[31] = rho; // b 




	/* ------------------------------------------*/
	/* --------------Flow restrictor-------------*/
	/* ------------------------------------------*/

	/*  Upstream parameters */
	double P0 = 0;
	double tLine = 304.724;
	double pLine = 0.8; // 0.6 Mpa;
	double pLine_FINAL = 3;
	double aReduced_final = 0.001;
	double molarMass = 29;
	double Rg = 8.314;
	double rLine = Rg/molarMass;
	double gammaLine = 1.4;
	double aReduced = 0.0003924;
	double vDead = (85*1000) ; /*  dead volume in mL -> mm^3 */

	/* ------------------------------------------*/
	/* ---------------Stretch rod----------------*/
	/* ------------------------------------------*/


	// stretch rod polynomial
	double a0 = -2.2264e7;
	double a1 = 2.3704e7;
	double a2 = -9.3769e6;
	double a3 = 1.6212e6;
	double a4 =-9.7380e4;
	double a5 = -1.8801e3;
	double a6 = 559.3131;
	double a7 = 0.2565;


	double stretchRodRad = 5.5;
	int numPointsRod = 15;

	MAT * srNodes = m_get(numPointsRod+2,2);

	MAT * srNodes_O = m_get(numPointsRod+2,2);
	for ( int i = 0 ; i < numPointsRod ; i++){
		double theta = -PI/2.00 + (PI/2/(numPointsRod-1))*i;
		srNodes->me[i][0] = stretchRodRad*cos(theta);
		// either 10.3 or 9
		srNodes->me[i][1] =10.1+stretchRodRad*sin(theta);
		srNodes_O->me[i][0] = srNodes->me[i][0];
		srNodes_O->me[i][1] = srNodes->me[i][1];
	}
	srNodes->me[numPointsRod][0] = stretchRodRad;
	srNodes->me[numPointsRod][1] = 70;
	srNodes_O->me[numPointsRod][0] = stretchRodRad;
	srNodes_O->me[numPointsRod][1] = 70;

	srNodes->me[numPointsRod+1][0] = 0;
	srNodes->me[numPointsRod+1][1] = 70;
	srNodes_O->me[numPointsRod+1][0] = 0;
	srNodes_O->me[numPointsRod+1][1] = 70;



	/* ------------------------------------------*/
	/* ---------------Mould----------------*/
	/* ------------------------------------------*/
	int numPoints_mould = 5;
	MAT * mould_Nodes = m_get(numPoints_mould,2);

	mould_Nodes->me[0][0] = 11.25;
	mould_Nodes->me[0][1] = 80;
	mould_Nodes->me[1][0] = 11.25;
	mould_Nodes->me[1][1] = 68.920;
	mould_Nodes->me[2][0] = 35.25;
	mould_Nodes->me[2][1] = 48.920;
	mould_Nodes->me[3][0] = mould_Nodes->me[2][0];
	mould_Nodes->me[3][1] = -105.7350;
	mould_Nodes->me[4][0] = 0;
	mould_Nodes->me[4][1] = mould_Nodes->me[3][1];



	/*////////////////////////////////////////////////////////// */
	/*////////////////////////////////////////////////////////// */


	/*////////////////////////////////////////////////////////// */
	/*                                                           */
	/*			        PREPRCOCESSING STAGE                     */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */


	// CREATE POINTS BASED ON DELAUNAY TRIANGULATION


	char opt[20] = "p";
	char fileName[30] = "preform";
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
	BOUNDING_BOX * bounding_box = create_bounding_box(-1, 40,
	-150, 100 , 0, 0);

	double cell_size[2] = {1.5,2};

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

	// for ( int i = 0 ; i < numnodes ; i++)
	// {
	// 	xI_copy->me[i][1] += -30;
	// 			xI_copy->me[i][0] += 0.2;

	// }

	// move_nodes(cells, &active_cells, cell_size, xI_copy);
	// move_nodes(cells, &active_cells, cell_size, xI_copy);


	// fp = fopen("nodes.csv","w");
	// for ( int i =0 ; i < numnodes ; i++)
	// {
	// 	for ( int k = 0 ; k < dim ; k++)
	// 	{
	// 		fprintf(fp, "%lf",xI_copy->me[i][k]);

	// 		if ( k < dim - 1)
	// 		{
	// 			fprintf(fp,",");
	// 		}

	// 	}
	// 	fprintf(fp, "\n");
	// }
	// fclose(fp);
	// write_active_cells("active_cells.csv",active_cells);

	// exit(0);
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
		.dmax_radial = dmax, .dmax_tensor = dmax_tensor,.temperatures = temperatures};
	
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
		write_material_points("materialpoints.csv", material_points);

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
	IVEC * neighbours;
	VEC * phi;

	VEC * inv_nodal_mass = v_get(mfree.num_nodes);

	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		inv_nodal_mass->ve[i] = 1.00/nodal_mass->ve[i];
	}

	printf("total mass = %lf \n", v_sum(nodal_mass));



	/* ------------------------------------------*/
	/* ------------Boundaries--------------------*/
	/* ------------------------------------------*/


	// Traction
	IVEC * traction_nodes ;
	getBoundary(&traction_nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,2);
	pressure_boundary * pB = new_pressure_boundary(traction_nodes, &mfree);
	m_foutput(stdout,pB->coords);
	pB->is_axi = is_AXI;

	// /*  EB1  */
	EBC * eb1 = malloc(1*sizeof(EBC));
	eb1->dofFixed = 3;
	getBoundary(&eb1->nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,5);
	//iv_addNode(eb1->nodes,traction_nodes->ive[0],'s');
	iv_addNode(eb1->nodes,traction_nodes->ive[traction_nodes->max_dim -1  ],'s');
	int num_nodes_eb1 = eb1->nodes->max_dim;
	setUpBC(eb1,inv_nodal_mass,&mfree);
	

	// /*  EB2 */
	EBC * eb2 = malloc(1*sizeof(EBC));
	eb2->dofFixed = 1;
	getBoundary(&eb2->nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,4);
	//iv_addNode(eb2->nodes,traction_nodes->ive[traction_nodes->max_dim - 1],'e');
	//iv_addNode(eb2->nodes,traction_nodes->ive[0],'e');

	int num_nodes_eb2 = eb2->nodes->max_dim;
	setUpBC(eb2,inv_nodal_mass,&mfree);

	m_foutput(stdout,eb1->coords);


	m_foutput(stdout,eb2->coords);


	// /*  EB3 */
	int numB3 = 13;
	IVEC * eb3_nodes = iv_get(numB3);
	MAT * contact_nodes_coords = m_get(numB3,dim);

	for ( int i = 0 ; i < numB3 ; i++){
		eb3_nodes->ive[i] = traction_nodes->ive[ i];
		contact_nodes_coords->me[i][0] = mfree.nodes->me[eb3_nodes->ive[i]][0];
		contact_nodes_coords->me[i][1] = mfree.nodes->me[eb3_nodes->ive[i]][1];

	}

		// get shape function and contact nodes
	shape_function_container * phi_contact = mls_shapefunction(contact_nodes_coords,1, &mfree);
	m_foutput(stdout,contact_nodes_coords);


	IVEC * eb4_nodes ;
	getBoundary(&eb4_nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,6);
	iv_addNode(eb4_nodes,eb2->nodes->ive[0],'e');
	int num_nodes_eb4 = eb4_nodes->max_dim;

	MAT * contact_mould_nodes_coords = m_get(num_nodes_eb4 ,dim);

	for ( int i = 0 ; i < num_nodes_eb4 ; i++){
		contact_mould_nodes_coords->me[i][0] = mfree.nodes->me[eb4_nodes->ive[i]][0];
		contact_mould_nodes_coords->me[i][1] = mfree.nodes->me[eb4_nodes->ive[i]][1];

	}
	// get shape function and contact nodes
	shape_function_container * phi_contact_mould = mls_shapefunction(contact_mould_nodes_coords,  1, &mfree);

	m_foutput(stdout,contact_mould_nodes_coords);




	/* ------------------------------------------*/
	/* -----------Transformation matrix----------*/
	/* ------------------------------------------*/


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


	/* ------------------------------------------*/
	/* --------------Cavity pressure-------------*/
	/* ------------------------------------------*/

	// /*  Cavity pressure paramters */
	double chokedMassRate =  flowRate(0.528, tLine, pLine*1e6, rLine, aReduced, gammaLine);
	double volumeInitial = cavityVolume(traction_nodes,mfree.nodes);
	double massInitial = (P0*volumeInitial)/(rLine*1000*tLine);
	double massAir = massInitial;
	double pRatio;



	// material
	double Kb = matParams->ve[9];
	double Gb = matParams->ve[10];


	double mu = Gb;
	double lambda = Kb - (2.00/3.00)*mu;


	// SPLIT MATERIAL POINTS
	int count;
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

	int total_number_material_points = 0;
	for (int i = 0 ; i < NUMBER_OF_THREADS ; i++)
	{
		total_number_material_points += material_point_groups[i]->max_dim;
	}
	printf("total number of material points = %d \n", total_number_material_points);

	for ( int i = 0 ; i < number_of_material_points ; i++)
	{

			material_points->MP[i]->stateNew->temperature=97.45+273.15;
			material_points->MP[i]->stateOld->temperature=97.45+273.15;
	
		

	}

	// v_zero(nodal_mass);
	// for ( int i = 0 ; i < number_of_material_points ; i++)
	// {
	// 	material_points->MP[i] = update_material_point(material_points->MP[i], cells,
	// 	mfree.nodes, nodal_mass);
	// }

	printf("mass = %lf \n", v_sum(nodal_mass));
	// ///////////////////////////////////////////////////////////////

	// /*////////////////////////////////////////////////////////// */
	// /*                                                           */
	// /*			      EXPLICIT TIME STEPPNG                      */
	// /*                                                           */
	// /*////////////////////////////////////////////////////////// */
	
	// time parameters
	double t_n = 0;
	double t_n_1 = 0;
	double t_n_h =  0; 

	int print_count = 0;



	int num_dof = dim*mfree.num_nodes;


	// External Forces
	VEC * Fext_n_1 = v_get(num_dof);
	VEC * Fext_n = v_get(num_dof);

	// Internal Forces
	VEC * Fint_n_1 = v_get(num_dof);
	VEC * Fint_n = v_get(num_dof);

	// Conctact Forces
	VEC * Fcont_n_1 = v_get(num_dof);


	// Net Force
	VEC * Fnet_n_1 = v_get(num_dof);

	// Kinematic variables
	// displacement
	VEC * d_n_1 = v_get(num_dof);
	VEC * d_n = v_get(num_dof);
	VEC * delta_disp = v_get(num_dof);
	VEC * nodal_disp = v_get(num_dof);
	// Velocity
	VEC * v_n_h = v_get(num_dof);
	VEC * v_n_mh = v_get(num_dof);
	VEC * v_n_1 = v_get(num_dof);
	VEC * v_n = v_get(num_dof);
	// Acceleration
	VEC * a_n_1 = v_get(num_dof);
	VEC * a_n = v_get(num_dof);
	VEC * D_N = v_get(num_dof);

	MAT * XI_n = m_get(mfree.num_nodes,2);
	MAT * XI_n_1 = m_get(mfree.num_nodes,2);
	VEC * R_pen = v_get(num_dof);

	m_copy(mfree.nodes,XI_n);

	/* Boundary conditions */
	eb1->uBar1 = v_get(eb1->nodes->max_dim);
	eb1->uBar2 = v_get(eb1->nodes->max_dim);
	eb2->uBar1 = v_get(eb2->nodes->max_dim);
	eb2->uBar1 = v_get(eb2->nodes->max_dim);


	VEC * v_correct = v_get(num_dof);
	// How often to write outputs
	int fileCounter = 0;

	// time step counter;
	int n = 1;

	// Energy
	double Wext = 0;
	double Wint = 0;
	double Wkin = 0;
	double Wbal = 0;

	// original position of the nodes
	MAT * nodes_X = m_copy(mfree.nodes,MNULL);

	// contact
	double f1Cor = 0;
	double f2Cor = 0; 
	MAT * msNormal = m_get(1,dim);
	MAT * testPoint = m_get(1,dim);
	double distanceProj = 0; 


	//External loading
	double pre_n = 0;
	double pre_n_1 = 0; 
	double volume = 0;
	double volume_t = 0;
	double pLine_n;

	double v_rod = 0;

	/*  For writing to file */
	fp = fopen("pressureTime.txt","w");
	fprintf(fp,"%lf %lf\n",0.0,0.0);
	fclose(fp);


	// stretch rod displacment 
	double disp_rod = 0;
	double disp_rod_n_1 = 0;
	double disp_rod_n = 0;

	// timing parameters
	struct timeval start3, end3;
	gettimeofday(&start3, NULL);


	// updated variables
	VEC * inc_disp = v_get(num_dof);


	VEC * FINT[NUMBER_OF_THREADS];
	VEC * RPEN[NUMBER_OF_THREADS];

	VEC * NODAL_MASS[NUMBER_OF_THREADS];
	for ( int i = 0 ; i < NUMBER_OF_THREADS ; i++)
	{

		NODAL_MASS[i] = v_get(numnodes);
		FINT[i] = v_get(num_dof);
		RPEN[i] = v_get(num_dof);

	}
	internal_force_args * INTERNAL_FORCE_ARGS = 
	calloc(NUMBER_OF_THREADS,sizeof(internal_force_args));

	// check if problem is axisymmetric
	int dim_piola;
	int dim_strain;
	int dim_cauchy;


	if ( is_AXI== 1){
		dim_piola = 5;
		dim_strain = 3;
		dim_cauchy = 4;
	}else{
		dim_piola = dim*dim;
		dim_strain = dim;
		dim_cauchy = (dim*dim) - (dim -1);
	}



	for ( int i = 0 ; i < NUMBER_OF_THREADS ; i++)
	{

		INTERNAL_FORCE_ARGS[i].NODAL_MASS = NODAL_MASS[i];
		INTERNAL_FORCE_ARGS[i].FINT = FINT[i];
		INTERNAL_FORCE_ARGS[i].RPEN = RPEN[i];
		INTERNAL_FORCE_ARGS[i].mat_points = material_point_groups[i];
		INTERNAL_FORCE_ARGS[i].material_points = material_points;
		INTERNAL_FORCE_ARGS[i].inc_disp = inc_disp;
		INTERNAL_FORCE_ARGS[i].materialParameters = matParams;
		INTERNAL_FORCE_ARGS[i].XI_n = XI_n;
		INTERNAL_FORCE_ARGS[i].XI_n_1 = XI_n_1;
		INTERNAL_FORCE_ARGS[i].cells = cells;
		INTERNAL_FORCE_ARGS[i].dt = delta_t;
		INTERNAL_FORCE_ARGS[i].G = m_get(dim_piola,dim_cauchy);
		INTERNAL_FORCE_ARGS[i].sigma = v_get(dim_cauchy);
		INTERNAL_FORCE_ARGS[i].velocity = v_n_h;
		INTERNAL_FORCE_ARGS[i].mass = nodal_mass;


	}



	/*  Explicit Loop */
	while ( t_n < TMAX)
	//while ( n < 10000)
	{

		/*  Update time step */
		t_n_1 = t_n + delta_t;


		for ( int k = 0 ; k < eb1->nodes->max_dim; k++)
		{
			int index = eb1->nodes->ive[k];
			a_n->ve[2*index] = 0;
			a_n->ve[2*index+1] = 0;

		}
		for ( int k = 0 ; k < eb2->nodes->max_dim ; k++)
		{
			int index = eb2->nodes->ive[k];
			a_n->ve[2*index] = 0;
		}



		/*  Make a time step  */ 
		__mltadd__(v_n_h->ve, a_n->ve,delta_t,num_dof);


		for ( int k = 0 ; k < eb1->nodes->max_dim; k++)
		{
			int index = eb1->nodes->ive[k];
			v_n_h->ve[2*index] = 0;
			v_n_h->ve[2*index+1] = 0;

		}
		for ( int k = 0 ; k < eb2->nodes->max_dim ; k++)
		{
			int index = eb2->nodes->ive[k];
			v_n_h->ve[2*index] = 0;
		}

		// for ( int k = 0 ; k < eb3_nodes->max_dim ; k++)
		// {
		// 	int index = eb3_nodes->ive[k];
		// 	v_n_h->ve[2*index+1] = -250*smoothstep(t_n_1,0.03,0);					

		// }


		__mltadd__(d_n_1->ve,v_n_h->ve,delta_t, num_dof);


		__add__(nodes_X->base, d_n_1->ve, XI_n_1->base, num_dof);


		/* ------------------------------------------*/
		/* -----------Contact Conditions-------------*/
		/* ------------------------------------------*/
		__zero__(Fcont_n_1->ve, num_dof);

#ifdef WITH_STRETCHROD
		__add__(nodes_X->base, d_n_1->ve, XI_n_1->base, num_dof);

		// IF STRETCH rod hasn't reached its maximum travel 
		if ( disp_rod_n < DISP_ROD_MAX){

		/*  Update stretch rod */
			double x = t_n_1*smoothstep(t_n_1,0.03,0);
			disp_rod_n_1 = a0*pow(x,7) + a1*pow(x,6) + a2*pow(x,5) + a3*pow(x,4) + a4*pow(x,3) + a5*pow(x,2) +a6*pow(x,1) + a7;
			disp_rod_n_1 = disp_rod_n_1;
			for ( int i = 0 ; i < srNodes->m ; i++){
	
				srNodes->me[i][1] = srNodes_O->me[i][1] - disp_rod_n_1;

			}
		}



		// STRETCH ROD CONTACT CONDITIONS

		for ( int i = 0 ; i < eb3_nodes->max_dim ; i++){

			neighbours = phi_contact->sf_list[i]->neighbours;
			phi = phi_contact->sf_list[i]->phi;

			int index = eb3_nodes->ive[i];
			testPoint->me[0][0] = XI_n_1->me[eb3_nodes->ive[i]][0];
			testPoint->me[0][1] = XI_n_1->me[eb3_nodes->ive[i]][1];
			distanceProj = contactDetection(testPoint,srNodes,msNormal);
			if (distanceProj > 0){


				f1Cor =1*(2*distanceProj*msNormal->me[0][0]*
					nodal_mass->ve[eb3_nodes->ive[i]])/pow(delta_t,2);
				f2Cor = 1*(2*distanceProj*msNormal->me[0][1]*
					nodal_mass->ve[eb3_nodes->ive[i]])/pow(delta_t,2);

				for ( int k = 0 ; k < neighbours->max_dim ; k++){
					Fcont_n_1->ve[2*neighbours->ive[k]] += phi->ve[k]*f1Cor; 
					Fcont_n_1->ve[2*neighbours->ive[k]+1] += phi->ve[k]*f2Cor; 
				}

			}
		}
#endif
#ifdef WITH_MOULD
		// MOULD CONTACT CONDITIONS
		for ( int i = 0 ; i < eb4_nodes->max_dim ; i++){


			neighbours = phi_contact_mould->sf_list[i]->neighbours;
			phi = phi_contact_mould->sf_list[i]->phi;
			testPoint->me[0][0] = XI_n_1->me[eb4_nodes->ive[i]][0];
			testPoint->me[0][1] = XI_n_1->me[eb4_nodes->ive[i]][1];

			distanceProj = contactDetection(testPoint,mould_Nodes,msNormal);

			if (distanceProj > 0){


				f1Cor = (2*distanceProj*msNormal->me[0][0]*nodal_mass->ve[eb4_nodes->ive[i]])/pow(delta_t,2);
				f2Cor = (2*distanceProj*msNormal->me[0][1]*nodal_mass->ve[eb4_nodes->ive[i]])/pow(delta_t,2);


				for ( int k = 0 ; k < neighbours->max_dim ; k++){
					Fcont_n_1->ve[2*neighbours->ive[k]] += phi->ve[k]*f1Cor; 
					Fcont_n_1->ve[2*neighbours->ive[k]+1] += phi->ve[k]*f2Cor; 
				}

			}


		}
#endif
#if defined(WITH_STRETCHROD) || defined (WITH_MOULD)
		// /*  Find a corrective acceleration - method in pronto 3D manual*/
		for ( int i = 0 ; i < numnodes  ; i++ )
		{
			a_n->ve[2*i] = Fcont_n_1->ve[2*i]/nodal_mass->ve[i];
			a_n->ve[2*i+1] = Fcont_n_1->ve[2*i+1]/nodal_mass->ve[i];
		}

		/* ------------------------------------------*/
		/* -----------Boundary Conditions------------*/
		/* ------------------------------------------*/




		__mltadd__(v_n_h->ve, a_n->ve,delta_t,num_dof);



		for ( int k = 0 ; k < eb2->nodes->max_dim ; k++)
		{
			int index = eb2->nodes->ive[k];
			v_n_h->ve[2*index] = 0;
		}



		__mltadd__(d_n_1->ve,v_n_h->ve,delta_t, num_dof);








#endif



		
		__add__(nodes_X->base, d_n_1->ve, XI_n_1->base, num_dof);

		#ifdef IS_UPDATED
			move_nodes(cells, &active_cells, cell_size, XI_n_1);
		#endif




		/* ------------------------------------------*/
		/* --------------Write outputs---------------*/
		/* ------------------------------------------*/

		// update nodal positions
		if (( n % WRITE_FREQ == 0 ) || (n == 1) ){

			char filename[50];
			// snprintf(filename, 50, "displacement_%d%s",fileCounter,".csv");
			// saveDisp(updatedNodes,state_n_1,"./Displacement",filename);
			snprintf(filename, 50, "displacement_%d%s",fileCounter,".txt");
			mat2csv(XI_n_1,"./Displacement",filename);

			snprintf(filename, 50, "srRod_%d%s",fileCounter,".csv");
			disp2csv(srNodes,"./srRod",filename);

			fp = fopen("pressureTime.txt","a");
			fprintf(fp,"%lf %lf\n",t_n_1,pre_n_1);
			fclose(fp);

			write_active_cells("active_cells.csv",active_cells);


			snprintf(filename, 50, "MaterialPoints/materialpoints_%d%s",fileCounter,".txt");
			write_material_points(filename, material_points);

			snprintf(filename, 50, "MaterialPoints/Domains/domains_%d%s",fileCounter,".txt");
			write_domains(filename, material_points);	



			/* ------------------------------------------*/
			/* --------------Print Outputs--------------*/
			/* ------------------------------------------*/
			state_variables * stateNew = material_points->MP[240]->stateNew;
			stateNew->F->me[2][1] = t_n_1;

			snprintf(filename, 50, "strain_%d%s",fileCounter,".txt");
			mat2csv(stateNew->F,"./History/Strain",filename);

			//snprintf(filename, 50, "Stress_%d%s",print_count,".txt");
			//mat2csv(stateNew[i]->sigma,"./History/Stress",filename);

			snprintf(filename, 50, "Bond_Stress_%d%s",fileCounter,".txt");
			mat2csv(stateNew->Sb,"./History/Stress",filename);
			snprintf(filename, 50, "Conformational_Stress_%d%s",fileCounter,".txt");
			mat2csv(stateNew->Sc,"./History/Stress",filename);
			stateNew->F->me[2][1] = 0;
			++print_count;



	








			fileCounter++;


		}






		/* ------------------------------------------*/
		/* ------------Find External Force-----------*/
		/* ------------------------------------------*/
		__zero__(Fext_n_1->ve,num_dof);
		/*  Find Cavity volume */
		volume = cavityVolume(traction_nodes,XI_n_1);


		if ( t_n_1 < 0.35)
		{
			pLine_n = pLine;
		}else{
			pLine_n = pLine + (pLine_FINAL-pLine) * smoothstep(t_n_1, 0.40, 0.350);
			aReduced = aReduced_final;
		}

		/*  Find Cavity pressure */
		pRatio = pre_n/pLine_n;
		if ( pRatio <= 0.528){
			massAir += chokedMassRate*delta_t;
		}else{
			massAir += flowRate(pRatio,tLine,pLine_n*1e6, rLine, aReduced,gammaLine)*delta_t;
		}
		pre_n_1 = ((P0*(volume - volumeInitial) + 1000*massAir*rLine*tLine)/(volume+vDead));

		/*  Update pressure load */
		update_pressure_boundary(pB, XI_n_1);
		assemble_pressure_load(Fext_n_1, pre_n_1, pB);







		/* ------------------------------------------*/
		/* ------------Find Internal Force-----------*/
		/* ------------------------------------------*/
		__zero__(R_pen->ve,num_dof);
		__zero__(Fint_n_1->ve,Fint_n_1->max_dim);


		// update incremental displacemnt
		__sub__(d_n_1->ve, D_N->ve, inc_disp->ve, num_dof);
		int i;


		for ( i = 0 ; i  < NUMBER_OF_THREADS ; i++)
		{
			__zero__(RPEN[i]->ve, num_dof);
			__zero__(FINT[i]->ve, num_dof);

		}

		#pragma omp parallel for num_threads(NUMBER_OF_THREADS)
		for (i=0; i < number_of_material_points; i++){

			int ID = omp_get_thread_num();

			INTERNAL_FORCE_ARGS[ID].MP = material_points->MP[i];
			internal_force_buckley(&INTERNAL_FORCE_ARGS[ID]);

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


		/* ------------------------------------------*/
		/* ---------------Find Net Force-------------*/
		/* ------------------------------------------*/

		/*  Balance of forces */
		__sub__(Fint_n_1->ve, R_pen->ve, Fint_n_1->ve, num_dof);
		__sub__(Fext_n_1->ve, Fint_n_1->ve, Fnet_n_1->ve,num_dof);



		/* -----------------------------------------------*/
		/* ---------------Update Acceleration-------------*/
		/* -----------------------------------------------*/

		for ( int i = 0 ; i < numnodes  ; i++ )
		{
			a_n_1->ve[2*i] = Fnet_n_1->ve[2*i]/nodal_mass->ve[i];
			a_n_1->ve[2*i+1] = Fnet_n_1->ve[2*i+1]/nodal_mass->ve[i];
		}
	
		/*  Integer time step velocity */
		v_mltadd(v_n_h,a_n_1,0.5*delta_t,v_n_1);	


		for ( int i = 0 ; i < num_nodes_eb1 ; i++)
		{

			v_n_1->ve[eb1->nodes->ive[i]*2] = 0;
			v_n_1->ve[eb1->nodes->ive[i]*2 + 1] = 0;
		}
		for ( int i = 0 ; i < num_nodes_eb2 ; i++)
		{

			v_n_1->ve[eb2->nodes->ive[i]*2] = 0;
		}


#ifdef IS_UPDATED
		if ( n % UPDATE_FREQUENCEY == 0){
			__zero__(nodal_mass->ve,numnodes);
			for ( i = 0 ; i  < NUMBER_OF_THREADS ; i++){
				__add__(NODAL_MASS[i]->ve, nodal_mass->ve,nodal_mass->ve, numnodes);
				__zero__(NODAL_MASS[i]->ve, numnodes);

			}

		}


#endif 


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

		/*  Update counters */
		t_n = t_n_1;
		/*  Store previous volume */
		volume_t = volume; 
		// Store previous time step quanities for the kinematic, and force variables.
		v_copy(Fint_n_1,Fint_n);
		v_copy(Fext_n_1,Fext_n);
		//v_copy(v_n_1,v_n_h);
		v_copy(d_n_1,d_n);
		//v_copy(v_n_1,v_n);
		v_copy(a_n_1,a_n);
		pre_n = pre_n_1;
		disp_rod_n = disp_rod_n_1;

		#ifdef IS_UPDATED	
			if ( n % UPDATE_FREQUENCEY == 0){
				m_copy(XI_n_1,XI_n);
				v_copy(d_n_1,D_N);
		}

		#endif		// update iteration counter
		double t_min = 0;
		n++	;

		if ( n % PRINT_FREQ == 0)
		{
			if ( print_count % 40 == 0)
			{
			printf("Iteration |    Time(s)  |   Energy    |  Pressure(MPa) | Volume (mL) | EXT  | MASS | INT |\n");

			}
			printf("%7i   |  %10.2E | %10.2E  |    %6lf    |   %7lf  |   %7lf  |  %7lf  |  %7lf    | \n",n,t_n,Wbal,
				pre_n_1,volume/1e3,v_sum(Fext_n_1),v_sum(nodal_mass),v_sum(Fint_n_1));
			++print_count;


		}





	}


	// /*////////////////////////////////////////////////////////// */
	// /*////////////////////////////////////////////////////////// */


	// gettimeofday(&endExp, NULL);

	// long elapsedExp = (endExp.tv_sec-endPre.tv_sec) + (endExp.tv_usec-endPre.tv_usec)/1000000.0;


	// printf("Post processing took %ld seconds \n Explicit routine took %ld seconds\n",elapsedPre,elapsedExp);




	// // /*////////////////////////////////////////////////////////// */
	// // /*			      POSTPROCESSING STAGE                       */
	// // /*                                                           */
	// // /*////////////////////////////////////////////////////////// */







	// /*////////////////////////////////////////////////////////// */
	// /*////////////////////////////////////////////////////////// */



	// /*  Program exit */
	exit(0);
}


static inline int heaviside(int x)
{

	if (x > 0){
		return 1;
	}else{
		return 0;
	}
}