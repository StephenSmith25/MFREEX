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

char * MATERIAL_NAME = "BUCKLEY";
const int BUCKLEY_MATERIAL = 1;
const int PLASTIC_MATERIAL = 0;

// time step parameters
const double TMAX = 1;
double delta_t = 5e-7;

// Meshfree parameters
const double dmax = 3;
const double dmax_x =2;
const double dmax_y =2;
double beta = 1.1;



char * basis_type = "linear";
char * weight = "quartic";
char * kernel_shape = "radial";


const int is_stabalised = 0;
const int is_constant_support_size = 1;

// stretch rod
const double DISP_ROD_MAX = 90; // 132;

// 
const int WITH_MOULD = 1;

const int WRITE_FREQ = 250;

int main(int argc, char** argv) {

	/*////////////////////////////////////////////////////////// */
	/*                                                           */
	/*			        Initialisation                           */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */

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
	matParams->ve[9] = 2.8e9; // Kb
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
		srNodes->me[i][1] =10.3+stretchRodRad*sin(theta);
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



	/*  Create points 
	 *
	 *	Via a 2D triangulation 
	 *
	 *  */
	// Read PLSG 

	char fileName[30] = "preform";
	double * points_out ;
	int * boundaryNodes = NULL;
	int numBoundary = 0;
	int * nodalMarkers;
	int  numnodes = 0;;
	double * temperatures;

	read_input_points(fileName,&points_out, &boundaryNodes, &nodalMarkers, &temperatures,&numnodes, &numBoundary);




	MAT * xI = m_get(numnodes,dim);
	// create nodes
	for ( int i = 0 ; i < numnodes ; i++)
	{
		xI->me[i][0] = points_out[2*i];
		xI->me[i][1] = points_out[2*i+1];

	}
	m_foutput(stdout, xI);
	fp = fopen("boundary.txt","w");
	for (int i = 0; i < numBoundary; ++i)
	{
		/* code */
		fprintf(fp,"%i\n",boundaryNodes[i]+1);
	}
	fclose(fp);


	
	struct timeval start, end;
	// generate clipped voronoi diagram
	gettimeofday(&start, NULL);
	voronoi_diagram * vor = NULL;
	vor = generate_voronoi(xI->base, boundaryNodes, xI->m, numBoundary, 2);
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



	/* ------------------------------------------*/
	/* ------------Meshfree Domain---------------*/
	/* ------------------------------------------*/

	// shape function parameters
	VEC * dI = v_get(xI->m);

	double dmax_tensor[dim];
	dmax_tensor[0] = dmax_x;
	dmax_tensor[1] = dmax_y;

	// meshfree domain
	meshfreeDomain mfree = {.nodes = xI, .di = dI, .num_nodes = xI->m, .dim = dim, .IS_AXI = is_AXI,
		.weight_function = weight, .kernel_shape = kernel_shape, 
		.basis_type = basis_type,.is_constant_support_size = is_constant_support_size,
		.dmax_radial = dmax, .dmax_tensor = dmax_tensor, .beta=beta};
	
	setDomain(&mfree);




	fp = fopen("domains.txt","w");
	if ( mfree.kernel_support == RADIAL)
	{
		for ( int i = 0 ; i < mfree.num_nodes ; i++)
			fprintf(fp,"%lf\n",mfree.di->ve[i]);

	}else if ( mfree.kernel_support == RECTANGULAR)
	{
		for ( int i = 0 ; i < mfree.num_nodes ; i++)
			fprintf(fp,"%lf,%lf\n",mfree.di_tensor->me[i][0],mfree.di_tensor->me[i][1]);
	}else if ( mfree.kernel_support == ELLIPTICAL)
	{		
			for ( int i = 0 ; i < mfree.num_nodes ; i++)
			fprintf(fp,"%lf,%lf,%lf,%lf\n",mfree.MI[i]->me[0][0],mfree.MI[i]->me[0][1],
				mfree.MI[i]->me[1][0],mfree.MI[i]->me[1][1]);

	}
	fclose(fp);
	
	/* ------------------------------------------*/
	/* ------------------SCNI--------------------*/
	/* ------------------------------------------*/

	SCNI_OBJ * _scni_obj = NULL;
	MSCNI_OBJ * _mscni_obj = NULL;

	struct timeval start2, end2;
	gettimeofday(&start2, NULL);
	// set up return container
	// generate shape functions at sample points 
	_scni_obj = generate_scni(vor, NULL , is_stabalised, is_AXI, dim, &mfree);
	gettimeofday(&end2, NULL);
	delta = ((end2.tv_sec  - start2.tv_sec) * 1000000u + 
         end2.tv_usec - start2.tv_usec) / 1.e6;
	printf("scni function took %lf seconds to run\n", delta);

	// Test divergence free condition
	IVEC * index ;
	MAT * B; 
	MAT * check_B = m_get(dim*dim,2);
	printf("checking scni \n \n");
	int checkPoint = 33;

	m_foutput(stdout,_scni_obj->scni[checkPoint]->B);
	iv_foutput(stdout, _scni_obj->scni[checkPoint]->sfIndex);

	printf("checking divergence free condition at point %d\n", checkPoint);
	printf("with coordinates %lf %lf\n", mfree.nodes->me[checkPoint][0], mfree.nodes->me[checkPoint][1]);
	printf("cell area = %lf \n", _scni_obj->scni[checkPoint]->area);
	printf("and neighbours \n");
	for ( int i = 0 ; i < _scni_obj->num_points ; i++)
	{
		index = _scni_obj->scni[i]->sfIndex;
		B = _scni_obj->scni[i]->B;

		for (int k = 0 ; k < index->max_dim ; k++){
			int indx = index->ive[k];

			if ( indx == checkPoint)
			{
				check_B->me[0][0] += B->me[0][2*k]*_scni_obj->scni[i]->area;
				check_B->me[1][1] += B->me[1][2*k+1]*_scni_obj->scni[i]->area;
				check_B->me[2][0] += B->me[2][2*k]*_scni_obj->scni[i]->area;
				check_B->me[3][1] += B->me[3][2*k+1]*_scni_obj->scni[i]->area;
			}
		}
	}
	// m_foutput(stdout, check_B);


	// m_foutput(stdout,_scni_obj->scni[checkPoint]->B);
	// iv_foutput(stdout, _scni_obj->scni[checkPoint]->sfIndex);
	


	/* ------------------------------------------*/
	/* ----------------Mass Vector---------------*/
	/* ------------------------------------------*/
	VEC * phi;
	IVEC * neighbours;
	VEC * nodal_mass = v_get(mfree.num_nodes);
	VEC * inv_nodal_mass = v_get(mfree.num_nodes);
	// get shape function and contact nodes
	shape_function_container * phi_nodes = mls_shapefunction(mfree.nodes, 1, &mfree);

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

		//nodal_mass->ve[i] = volume*rho;
		//inv_nodal_mass->ve[i] = 1.000/nodal_mass->ve[i];

	}
	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		inv_nodal_mass->ve[i] = 1.00/nodal_mass->ve[i];
	}
	printf("total mass = %lf (g) \n", v_sum(nodal_mass)*1000);
	
	/* ------------------------------------------*/
	/* ------------Boundaries--------------------*/
	/* ------------------------------------------*/

	m_foutput(stdout, xI);

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
	//iv_addNode(eb1->nodes,traction_nodes->ive[traction_nodes->max_dim -1  ],'s');

	int num_nodes_eb1 = eb1->nodes->max_dim;
	setUpBC(eb1,inv_nodal_mass,&mfree);
	

	// /*  EB2 */
	EBC * eb2 = malloc(1*sizeof(EBC));
	eb2->dofFixed = 1;
	getBoundary(&eb2->nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,4);
	//iv_addNode(eb2->nodes,traction_nodes->ive[traction_nodes->max_dim - 1],'e');
	iv_addNode(eb2->nodes,traction_nodes->ive[0],'e');

	int num_nodes_eb2 = eb2->nodes->max_dim;
	setUpBC(eb2,inv_nodal_mass,&mfree);

	m_foutput(stdout,eb1->coords);


	m_foutput(stdout,eb2->coords);



	// /*  EB3 */
	int numB3 = 15;
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

	/* ------------------------------------------*/
	/* --------------State storage---------------*/
	/* ------------------------------------------*/
	state_variables ** state_n = new_material_states(temperatures, mfree.num_nodes, 
		BUCKLEY_MATERIAL , 
		PLASTIC_MATERIAL , 
		dim, is_AXI);
	state_variables ** state_n_1 = new_material_states(temperatures, mfree.num_nodes, 
		BUCKLEY_MATERIAL , 
		PLASTIC_MATERIAL , 
		dim, is_AXI);
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
	MAT * updatedNodes = m_get(mfree.num_nodes,2);


	/* Boundary conditions */
	eb1->uBar1 = v_get(eb1->nodes->max_dim);
	eb1->uBar2 = v_get(eb1->nodes->max_dim);
	eb2->uBar1 = v_get(eb2->nodes->max_dim);
	eb2->uBar1 = v_get(eb2->nodes->max_dim);


	VEC * v_correct = v_get(num_dof);
	// How often to write outputs
	int fileCounter = 0;

	// time step counter;
	int n = 0;

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

	VEC * disp_inc = v_get(num_dof);
	VEC * disp_r = v_get(num_dof);



	/*  Explicit Loop */
	while ( t_n < TMAX)
	//while ( n < 1)
	{

		/*  Update time step */
		t_n_1 = t_n + delta_t;
		/*  Make a time step  */ 
		__mltadd__(v_n_h->ve, a_n->ve,delta_t,num_dof);
		__mltadd__(d_n_1->ve,v_n_h->ve,delta_t, num_dof);


		mv_mlt(Lambda,d_n_1,nodal_disp);
		__add__(nodes_X->base, nodal_disp->ve, updatedNodes->base, num_dof);

		/* ------------------------------------------*/
		/* -----------Contact Conditions-------------*/
		/* ------------------------------------------*/
		__zero__(Fcont_n_1->ve, num_dof);

		// STRETCH ROD 
		// update stretch rod position
		if ( disp_rod_n < DISP_ROD_MAX){
		/*  Update stretch rod */
			double x = t_n_1*smoothstep(t_n_1,0.001,0);

			disp_rod_n_1 = a0*pow(x,7) + a1*pow(x,6) + a2*pow(x,5) + a3*pow(x,4) + a4*pow(x,3) + a5*pow(x,2) +a6*pow(x,1) + a7;
			disp_rod_n_1 = disp_rod_n_1;
	
			for ( int i = 0 ; i < srNodes->m ; i++){
	
				//srNodes->me[i][1] = srNodes_O->me[i][1] + disp_rod_n_1;
				srNodes->me[i][1] = srNodes_O->me[i][1] - disp_rod_n_1;

			}
		}
		// STRETCH ROD CONTACT CONDITIONS

		for ( int i = 0 ; i < eb3_nodes->max_dim ; i++){

			neighbours = phi_contact->sf_list[i]->neighbours;
			phi = phi_contact->sf_list[i]->phi;
			testPoint->me[0][0] = updatedNodes->me[eb3_nodes->ive[i]][0];
			testPoint->me[0][1] = updatedNodes->me[eb3_nodes->ive[i]][1];
			distanceProj = contactDetection(testPoint,srNodes,msNormal);
			if (distanceProj > 0){

				f1Cor = (2*distanceProj*msNormal->me[0][0]*nodal_mass->ve[eb3_nodes->ive[i]])/pow(delta_t,2);
				f2Cor = (2*distanceProj*msNormal->me[0][1]*nodal_mass->ve[eb3_nodes->ive[i]])/pow(delta_t,2);


				for ( int k = 0 ; k < neighbours->max_dim ; k++){
					Fcont_n_1->ve[2*neighbours->ive[k]] += phi->ve[k]*f1Cor; 
					Fcont_n_1->ve[2*neighbours->ive[k]+1] += phi->ve[k]*f2Cor; 
				}

			}


		}
		if ( WITH_MOULD == 1){
		// MOULD CONTACT CONDITIONS
		for ( int i = 0 ; i < eb4_nodes->max_dim ; i++){

			neighbours = phi_contact_mould->sf_list[i]->neighbours;
			phi = phi_contact_mould->sf_list[i]->phi;
			testPoint->me[0][0] = updatedNodes->me[eb4_nodes->ive[i]][0];
			testPoint->me[0][1] = updatedNodes->me[eb4_nodes->ive[i]][1];

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
		}


		// /*  Find a corrective acceleration - method in pronto 3D manual*/
		for ( int i = 0 ; i < numnodes  ; i++ )
		{
			a_n->ve[2*i] = Fcont_n_1->ve[2*i]*inv_nodal_mass->ve[i];
			a_n->ve[2*i+1] = Fcont_n_1->ve[2*i+1]*inv_nodal_mass->ve[i];
		}

		__mltadd__(v_n_h->ve, a_n->ve,delta_t,num_dof);
		__mltadd__(d_n_1->ve,v_n_h->ve,delta_t, num_dof);


		/* ------------------------------------------*/
		/* -----------Boundary Conditions------------*/
		/* ------------------------------------------*/
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

		// Symmetry boundary /
		enforceBC(eb2,d_n_1); 
		sv_mlt(1.00/(delta_t),eb2->uCorrect1,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k] += v_correct->ve[k];
		}

		// find new nodal positions
		mv_mlt(Lambda,d_n_1,nodal_disp);
		__add__(nodes_X->base, nodal_disp->ve, updatedNodes->base, num_dof);

	
		/* ------------------------------------------*/
		/* ------------Find External Force-----------*/
		/* ------------------------------------------*/
		__zero__(Fext_n_1->ve,num_dof);
		/*  Find Cavity volume */
		volume = cavityVolume(traction_nodes,updatedNodes);


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
		update_pressure_boundary(pB, updatedNodes);

		assemble_pressure_load(Fext_n_1, pre_n_1, pB);


		/* ------------------------------------------*/
		/* ------------Find Internal Force-----------*/
		/* ------------------------------------------*/
		__sub__(d_n_1->ve, disp_r->ve,disp_inc->ve, num_dof);

		internalForce_Inelastic_Buckley(Fint_n_1, _scni_obj,
		disp_inc, v_n_h,
		matParams, state_n_1, state_n,
		mfree.IS_AXI, dim,delta_t,t_n_1, MATERIAL_NAME);


		/* ------------------------------------------*/
		/* ---------------Find Net Force-------------*/
		/* ------------------------------------------*/

		/*  Balance of forces */
		__sub__(Fext_n_1->ve, Fint_n_1->ve, Fnet_n_1->ve,num_dof);


		/* -----------------------------------------------*/
		/* ---------------Update Acceleration-------------*/
		/* -----------------------------------------------*/

		for ( int i = 0 ; i < numnodes  ; i++ )
		{
			a_n_1->ve[2*i] = Fnet_n_1->ve[2*i]*inv_nodal_mass->ve[i];
			a_n_1->ve[2*i+1] = Fnet_n_1->ve[2*i+1]*inv_nodal_mass->ve[i];
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


		/* ------------------------------------------*/
		/* --------------Write outputs---------------*/
		/* ------------------------------------------*/

		// update nodal positions
		if ( n % WRITE_FREQ == 0 ){

			char filename[50];
			snprintf(filename, 50, "displacement_%d%s",fileCounter,".csv");
			saveDisp(updatedNodes,state_n_1,"./Displacement",filename);

			snprintf(filename, 50, "srRod_%d%s",fileCounter,".csv");
			disp2csv(srNodes,"./srRod",filename);

			fp = fopen("pressureTime.txt","a");
			fprintf(fp,"%lf %lf\n",t_n_1,pre_n_1);
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
		// update iteration counter
		double t_min = 0;
		n++	;



		if ( n % WRITE_FREQ == 0)
			printf("%i  \t  %lf %10.2E %lf %lf %lf %lf %10.2E \n",n,t_n,Wbal,
				pre_n_1,disp_rod_n,v_rod,volume/1e3,delta_t);




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