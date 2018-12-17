#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "matrix2.h"
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
#include "Force/Internal/internalForce_hyperelastic_S.h"
#include "mat2csv.h"
#include "trigen.h"
#include "Boundary/getBoundary.h"
#include "Boundary/iv_addNode.h"
#include "Boundary/Traction/Pressure/new_pressure_load.h"
#include "Boundary/Traction/Cavity/cavityVolume.h"
#include "Boundary/Traction/Cavity/flowRate.h"
#include <math.h>
#include "Deformation/poldec.h"
#include "Material/Buckley/new_Buckley_State.h"


// #include "internalForceBuckley.h"
// #include "contactDetection.h"



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

	/*  Material struct */
	VEC * matParams = v_get(26);
	matParams->ve[0] = 2.814e-3; // VS
	matParams->ve[1] = 0.526e-3; // VP
	matParams->ve[2] = (1.71e6); // mu*_0
	matParams->ve[3] = (328.76); // Tinf
	matParams->ve[4] = 358.15; // T*
	matParams->ve[5] = matParams->ve[4]; // Tf*
	matParams->ve[6] = (67.47); // Cv
	matParams->ve[7] = 1.23e5; // H0
	matParams->ve[8] = 8.314; // R
	matParams->ve[9] = 1.8e9; // Kb
	matParams->ve[10] = 6e8;// Gb
	matParams->ve[11] = -100000000000;// temperature
	matParams->ve[12] = matParams->ve[11];// Tf
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
	VEC * critLambdaParams = v_get(6);
	critLambdaParams->ve[0] = -0.0111; // C1
	critLambdaParams->ve[1] = 3.627; // C2
	critLambdaParams->ve[2] = 0.9856; // BETA
	critLambdaParams->ve[3] = -0.0356; // k
	critLambdaParams->ve[4] = 15.393; // b 
	critLambdaParams->ve[5] = matParams->ve[11]; // temperature 

	/*  Time step variables */

	double deltaT = 2e-7;
	double tMax = 0.8;





	/* ------------------------------------------*/
	/* --------------Flow restrictor-------------*/
	/* ------------------------------------------*/

	/*  Upstream parameters */
	double P0 = 0;
	double tLine = 304.724;
	double pLine = 0.8; // 0.6 Mpa;
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

	MAT * srNodes = m_get(numPointsRod+1,2);

	MAT * srNodes_O = m_get(numPointsRod+1,2);
	for ( int i = 0 ; i < numPointsRod ; i++){
		double theta = -PI/2.00 + (PI/2/(numPointsRod-1))*i;
		srNodes->me[i][0] = stretchRodRad*cos(theta);
		srNodes->me[i][1] = 9.21+stretchRodRad*sin(theta);
		srNodes_O->me[i][0] = srNodes->me[i][0];
		srNodes_O->me[i][1] = srNodes->me[i][1];
	}
	srNodes->me[numPointsRod][0] = stretchRodRad;
	srNodes->me[numPointsRod][1] = 70;
	srNodes_O->me[numPointsRod][0] = stretchRodRad;
	srNodes_O->me[numPointsRod][1] = 70;

	/*  Stretch rod */
	double dRodStop = 135; /*  130 mm  */
	double tStop = 0.5;
	double tRampRod = 1e-9;
	double vRod = 0;

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

	char opt[20] = "pDq0a2";
	char fileName[30] = "preform";
	double * points_out ;
	int * boundaryNodes;
	int numBoundary;
	int * nodalMarkers;
	int  numnodes;
	double * temperatures;
	trigen(&points_out,&boundaryNodes,opt,fileName,&numnodes,&numBoundary,&nodalMarkers,&temperatures);	


	MAT * xI = m_get(numnodes,dim);
	// create nodes
	for ( int i = 0 ; i < numnodes ; i++)
	{
		xI->me[i][0] = points_out[2*i];
		xI->me[i][1] = points_out[2*i+1];

	}
	
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
	double dmax = 2;
	int constant_support_size = 1;
	VEC * dI = v_get(xI->m);

	// meshfree domain
	meshfreeDomain mfree = {.nodes = xI, .di = dI, .num_nodes = xI->m, .dim = dim};
	setDomain(&mfree,constant_support_size, dmax);


	/* ------------------------------------------*/
	/* ------------Boundaries--------------------*/
	/* ------------------------------------------*/


	// Traction
	IVEC * traction_nodes ;
	getBoundary(&traction_nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,2);
	pressure_boundary * pB = new_pressure_load(traction_nodes, &mfree);


	// /*  EB1  */
	IVEC * eb1_nodes;
	getBoundary(&eb1_nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,5);
	iv_addNode(eb1_nodes,traction_nodes->ive[0],'s');

	// /*  EB2 */
	IVEC * eb2_nodes;
	getBoundary(&eb2_nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,4);
	iv_addNode(eb2_nodes,traction_nodes->ive[traction_nodes->max_dim - 1],'e');


	// /*  EB3 */
	int numB3 = 15;
	IVEC * eb3_nodes = iv_get(numB3);
	for ( int i = 0 ; i < numB3 ; i++){
		eb3_nodes->ive[i] = traction_nodes->ive[traction_nodes->max_dim -1 - i];
	}


	// print coordinates of boundary nodes
	// m_foutput(stdout,pB->coords);


	// printf("coordinates of eb1 nodes = \n");


	// for ( int i = 0 ; i < eb1_nodes->max_dim ; i++)
	// {
	// 	printf("%lf %lf \n", mfree.nodes->me[eb1_nodes->ive[i]][0], mfree.nodes->me[eb1_nodes->ive[i]][1]);
	// }


	// printf("coordinates of eb2 nodes = \n");

	// for ( int i = 0 ; i < eb2_nodes->max_dim ; i++)
	// {
	// 	printf("%lf %lf \n", mfree.nodes->me[eb2_nodes->ive[i]][0], mfree.nodes->me[eb2_nodes->ive[i]][1]);
	// }
	// printf("coordinates of eb3 nodes = \n");

	// for ( int i = 0 ; i < eb3_nodes->max_dim ; i++)
	// {
	// 	printf("%lf %lf \n", mfree.nodes->me[eb3_nodes->ive[i]][0], mfree.nodes->me[eb3_nodes->ive[i]][1]);
	// }
	// MAT * point = m_get(1,2);
	// EFG_SF * phiContact = malloc(eb3->nodes->max_dim * sizeof(EFG_SF));

	// define_support(phiContact,eb3->coords,efgBlock);
	// MLS_shapefunction(eb3->coords,phiContact,efgBlock);


	/* ------------------------------------------*/
	/* ------------------SCNI--------------------*/
	/* ------------------------------------------*/

	int is_stabalised = 0;
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
	m_foutput(stdout, check_B);
	


	/* ------------------------------------------*/
	/* ----------------Mass Vector---------------*/
	/* ------------------------------------------*/

	VEC * nodal_mass = v_get(mfree.num_nodes);
	VEC * inv_nodal_mass = v_get(mfree.num_nodes);

	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		double volume = _scni_obj->scni[i]->area;
		if ( is_AXI == 1)
		{
			double r = _scni_obj->scni[i]->r;
			volume = volume*r*2*PI;
		}
		nodal_mass->ve[i] = volume*rho;
		inv_nodal_mass->ve[i] = 1.000/nodal_mass->ve[i];

	}

	

	/* ------------------------------------------*/
	/* -----------Transformation matrix----------*/
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

	/* ------------------------------------------*/
	/* --------------Cavity pressure-------------*/
	/* ------------------------------------------*/

	// /*  Cavity pressure paramters */
	double chokedMassRate =  flowRate(0.528, tLine, pLine*1e6, rLine, aReduced, gammaLine);
	double volumeInitial = cavityVolume(traction_nodes,mfree.nodes);
	double massInitial = (P0*volumeInitial)/(rLine*1000*tLine);
	double massAir = massInitial;
	double pRatio;


	// tst poldec

	int test_dim = 3;
	MAT * F = m_get(test_dim,test_dim);
	F->me[0][0] = 1;
	F->me[0][1] = 0.495;
	F->me[0][2] = 0.5;

	F->me[1][0] = -0.333;
	F->me[1][1] = 1;
	F->me[1][2] = -0.247;


	F->me[2][0] = 0.959;
	F->me[2][1] = 0;
	F->me[2][2] = 1.5;


	MAT * R = m_get(test_dim,test_dim);
	MAT * U = m_get(test_dim,test_dim);
	MAT * V = m_get(test_dim,test_dim);


	poldec(F, R, U, V);
	printf("F = \n");
	m_foutput(stdout,F);
	printf("R = \n");
	m_foutput(stdout,R);
	printf("U = \n");
	m_foutput(stdout,U);
	printf("V = \n");
	m_foutput(stdout,V);

	/* ------------------------------------------*/
	/* --------------State storage---------------*/
	/* ------------------------------------------*/
	// store F, D at each time step for each material point for use with buckley model 
	state_Buckley ** state_n = new_Buckley_State(mfree.num_nodes,NULL,is_AXI,dim);;
	state_Buckley ** state_n_1 = new_Buckley_State(mfree.num_nodes,NULL,is_AXI,dim);


	// ///////////////////////////////////////////////////////////////

	// /*////////////////////////////////////////////////////////// */
	// /*                                                           */
	// /*			      EXPLICIT TIME STEPPNG                      */
	// /*                                                           */
	// /*////////////////////////////////////////////////////////// */
	
	// time parameters
	double t_max = 0.5; // 1s
	double delta_t = 1.55e-6;
	double t_n = 0;
	double t_n_1 = 0;


	int num_dof = dim*mfree.num_nodes;


	// External Forces
	VEC * Fext_n_1 = v_get(num_dof);
	VEC * Fext_n = v_get(num_dof);

	// Internal Forces
	VEC * Fint_n_1 = v_get(num_dof);
	VEC * Fint_n = v_get(num_dof);

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

	// original position of the nodes
	MAT * nodes_X = m_copy(mfree.nodes,MNULL);


	//External loading
	double pre_n = 0;
	double pre_n_1 = 0; 
	double volume = 0;
	double volume_t = 0;

	/*  For writing to file */
	fp = fopen("pressureTime.txt","w");
	fprintf(fp,"%lf %lf\n",0.0,0.0);
	fclose(fp);


	// /*  Contact with stretch rod */
	// double dispRod = 0;
	// double dispRodStop = 0;
	// double distanceProj = 0;
	// MAT * msNormal = m_get(1,2);
	// double f1Cor = 0;
	// double f2Cor = 0;

	// double beta_1 = 10;
	// double beta_2 = 10;



	// mv_mlt(Lambda,d_n_1,nodalDisp);
	// for ( int i = 0 ; i < efgBlock->numnode ; i++){
	// 	updatedNodes->me[i][0]  = efgBlock->nodes->me[i][0] + nodalDisp->ve[2*i] ;
	// 	updatedNodes->me[i][1]  = efgBlock->nodes->me[i][1] + nodalDisp->ve[2*i+1] ;	

	// 	updatedNodesT->me[i][0]  = efgBlock->nodes->me[i][0] + nodalDisp->ve[2*i] ;
	// 	updatedNodesT->me[i][1]  = efgBlock->nodes->me[i][1] + nodalDisp->ve[2*i+1] ;	

	// }

	struct timeval start3, end3;
	gettimeofday(&start3, NULL);



	// double disp_rod = 0;
	// /*  Explicit Loop */
	// while ( t_n < tMax){

	// 	/*  Update time step */
	// 	t_n_1 = t_n + deltaT;
	// 	t_n_h = (1.00/2.00)*(t_n + t_n_1);

	// 	/*  Make a time step  */ 
	// 	v_mltadd(v_n,a_n,(double)0.5*deltaT,v_n_h);
	// 	v_mltadd(d_n,v_n_h,deltaT,d_n_1);

	// 	if ( disp_rod < 132){
	// 	/*  Update stretch rod */
	// 		disp_rod = a0*pow(t_n_1,7) + a1*pow(t_n_1,6) + a2*pow(t_n_1,5) + a3*pow(t_n_1,4) + a4*pow(t_n_1,3) + a5*pow(t_n_1,2) +
	// 		a6*pow(t_n_1,1) + a7;
	// 		//disp_rod = vRod*t_n_1;

	// 		for ( int i = 0 ; i < srNodes->m ; i++){
	
	// 			srNodes->me[i][1] = srNodes_O->me[i][1] - disp_rod;

	// 		}
	// 		dispRod = disp_rod;
	// 	}
	// 	/*  Stretch rod boundary conditions */

	// 	for ( int i = 0 ; i < eb3->nodes->max_dim ; i++){


	// 		testPoint->me[0][0] = updatedNodes->me[eb3->nodes->ive[i]][0];
	// 		testPoint->me[0][1] = updatedNodes->me[eb3->nodes->ive[i]][1];

	// 		distanceProj = contactDetection(testPoint,srNodes,msNormal);

	// 		//printf("distance = %lf with normal {%lf,%lf}\n",distanceProj,msNormal->me[0][0],msNormal->me[0][1]);


	// 		if (distanceProj > 0){

	// 			f1Cor = 1*(2*distanceProj*msNormal->me[0][0]*mass->me[eb3->nodes->ive[i]*2][eb3->nodes->ive[i]*2])/pow(deltaT,2);
	// 			f2Cor = 1*(2*distanceProj*msNormal->me[0][1]*mass->me[eb3->nodes->ive[i]*2][eb3->nodes->ive[i]*2])/pow(deltaT,2);



	// 			for ( int k = 0 ; k < phiContact[i].index->max_dim ; k++){
	// 				Fnet->ve[2*phiContact[i].index->ive[k]] += phiContact[i].phi->ve[k]*f1Cor; 
	// 				Fnet->ve[2*phiContact[i].index->ive[k]+1] += phiContact[i].phi->ve[k]*f2Cor; 
	// 			}

	// 		}


	// 			// penalty contact
	// 			//double gamma_N = vRod - v_n_h->ve[eb3->nodes->ive[i]*2];
	// 			//double g_N = distanceProj;
	// 			// penalty method
	// 			//double p_bar = g_N gamma_N*beta_2*heaviside(gamma_N)


	// 	}
	// 	mv_mlt(invMass,Fnet,a_n);
	// 	/*  Find a corrective acceleration - method in pronto 3D manual*/
	// 	/*  Make a time step  */ 
	// 	v_mltadd(v_n,a_n,(double)0.5*deltaT,v_n_h);
	// 	v_mltadd(d_n,v_n_h,deltaT,d_n_1);


	// 	/*  Implement BCs */
	// 	enforceBC(eb1,d_n_1); 
	// 	// find velocity correction
	// 	sv_mlt(1.00/(2*deltaT),eb1->uCorrect1,v_correct);
	// 	for ( int k = 0 ; k < v_correct->max_dim; k++){
	// 		v_n_h->ve[2*k] += v_correct->ve[k];
	// 	}

	// 	sv_mlt(1.000/(2*deltaT),eb1->uCorrect2,v_correct);
	// 	for ( int k = 0 ; k < v_correct->max_dim; k++){
	// 		v_n_h->ve[2*k+1] += v_correct->ve[k];
	// 	}


	// 	// Symmetry boundary /
	// 	enforceBC(eb2,d_n_1); 
	// 	sv_mlt(1.00/(2*deltaT),eb2->uCorrect1,v_correct);
	// 	for ( int k = 0 ; k < v_correct->max_dim; k++){
	// 		v_n_h->ve[2*k] += v_correct->ve[k];
	// 	}


	// 	/*  Update nodal poistions */
	// 	mv_mlt(Lambda,d_n_1,nodalDisp);
	// 	for ( int i = 0 ; i < efgBlock->numnode ; i++){
	// 		updatedNodes->me[i][0]  = efgBlock->nodes->me[i][0] + nodalDisp->ve[2*i] ;
	// 		updatedNodes->me[i][1]  = efgBlock->nodes->me[i][1] + nodalDisp->ve[2*i+1] ;

	// 		if ( n % 10 == 0 )
	// 			{
	// 			updatedNodesT->me[i][0]  = efgBlock->nodes->me[i][0] + nodalDisp->ve[2*i] ;
	// 			updatedNodesT->me[i][1]  = efgBlock->nodes->me[i][1] + nodalDisp->ve[2*i+1] ;	
	// 		}

	// 	}

	// 	// if ( n % 100000 == 0) 
	// 	// {	
	// 	// 	 gpc_polygon ** voronoi_1;
	// 	// 	 callVoronoi(&voronoi_1,updatedNodes->base,numnodes);
 // 	// 		 clipVoronoi(&voronoi_1,updatedNodes->base,boundaryNodes,numnodes,numBoundary);
	// 	// }


	// 	/*  Find Cavity volume */
	// 	volume = cavityVolume(traction->nodes,updatedNodes);

	// 	if ( t_n_1 < t_pressure)
	// 	{
	// 		pre_n_1 = 0;
	// 	}else{


	// 	/*  Find Cavity pressure */
	// 	pRatio = pre_n/pLine;
	// 	if ( pRatio <= 0.528){
	// 		massAir += chokedMassRate*deltaT;
	// 	}else{
	// 		massAir += flowRate(pRatio,tLine,pLine*1e6, rLine, aReduced,gammaLine)*deltaT;
	// 	}
	// 	pre_n_1 = ((P0*(volume - volumeInitial) + 1000*massAir*rLine*tLine)/(volume+vDead));
	// 	}


	// 	//pre_n_1 = 0.5*smoothstep(t_n_1,0.2,t_pressure);

	// 	/*  Update pressure load */
	// 	traction->pressure = -pre_n_1;
	// 	externalForce(Fext_n_1,traction,scni,updatedNodes);


	// 	/*  Internal force */
	// 	//internalForce(Fint_n_1,&scni,d_n_1,matParams,efgBlock->numnode);
	// 	internalForceBuckley(Fint_n_1,scni,d_n_1,matParams,critLambdaParams,state_n,deltaT,efgBlock->numnode,t_n_1);

	// 	/*  Damping */
	// 	mv_mlt(mass,v_n_h,Fdamp);
	// 	sv_mlt(alpha,Fdamp,Fdamp);


	// 	/*  Balance of forces */
	// 	v_sub(Fext_n_1,Fint_n_1,Fnet);
	// 	v_sub(Fnet,Fdamp,Fnet);

	// 	/*  Find acceleration */
	// 	mv_mlt(invMass,Fnet,a_n_1);

	// 	/*  Integer time step velocity */
	// 	v_mltadd(v_n_h,a_n_1,t_n_1-t_n_h,v_n_1);	


	// 	// update nodal positions
	// 	if ( n % writeFreq == 0 ){
	// 		char * filename[50];
	// 		snprintf(filename, 50, "displacement_%d%s",fileCounter,".txt");
	// 		mat2csv(updatedNodes,"./Displacement",filename);

	// 		snprintf(filename, 50, "srRod_%d%s",fileCounter,".txt");
	// 		mat2csv(srNodes,"./srRod",filename);

	// 		fp = fopen("pressureTime.txt","a");
	// 		fprintf(fp,"%lf %lf\n",t_n_1,pre_n_1);
	// 		fclose(fp);

	// 		fileCounter++;
	// 	}


	// 	/*  Find energy */
	// 	v_sub(d_n_1,d_n,deltaDisp);
	// 	// External Energy
	// 	Wext_n_1 = Wext_n + 0.5*(in_prod(deltaDisp,Fext_n) + in_prod(deltaDisp,Fext_n_1));
	// 	// Internl Energy
	// 	Wint_n_1 = Wint_n + 0.5*(in_prod(deltaDisp,Fint_n) + in_prod(deltaDisp,Fint_n_1));
	// 	// Kinetic Energy = 1/2 v^T * M * v
	// 	Wkin_n_1 = 0 ;
	// 	for (int i = 0 ; i < 2* efgBlock->numnode ; i++){

	// 		Wkin_n_1 = Wkin_n_1 + 0.5*(v_n_1->ve[i]*mass->me[i][i]*v_n_1->ve[i]);
	// 	}
	// 	// Find the energy balance.
	// 	Wbal = fabs(Wkin_n_1 + Wint_n_1 - Wext_n_1) ;
	// 	// Find maximum of the energies, used to normalise the balance. Uses the max macro defined in matrix.h
	// 	double W_max = max(Wext_n_1,Wint_n_1);
	// 	W_max = max(W_max,Wkin_n_1);
	// 	/*  Update counters */
	// 	t_n = t_n_1;
	// 	// Store energy accumulated
	// 	Wext_n = Wext_n_1;
	// 	Wint_n = Wint_n_1;
	// 	/*  Store previous volume */
	// 	volume_t = volume; 
	// 	// Store previous time step quanities for the kinematic, and force variables.
	// 	v_copy(Fint_n_1,Fint_n);
	// 	v_copy(Fext_n_1,Fext_n);
	// 	v_copy(v_n_h,v_n_mh);
	// 	v_copy(d_n_1,d_n);
	// 	v_copy(v_n_1,v_n);
	// 	v_copy(a_n_1,a_n);
	// 	pre_n = pre_n_1;
	// 	// update iteration counter
	// 	n++	;
	// 	printf("%i  \t  %lf %10.2E %lf %lf %lf \n",n,t_n,Wbal/W_max,pre_n_1,dispRod,volume/1e3);




	// }


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
	return 0;
}


static inline int heaviside(int x)
{

	if (x > 0){
		return 1;
	}else{
		return 0;
	}
}