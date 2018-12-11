#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structure.h"
#include "vGeom.h"
#include "triangle.h"
#include "setDomain.h"
#include "matrix.h"
#include "matrix2.h"
#include "generateSCNI.h"
#include "getBoundary.h"
#include <trigen.h>
#include "iv_addNode.h"
#include "setUpBC.h"
#include "jc_voronoi.h"
#include "clipVoronoi.h"
#include "externalForce.h"
#include <pthread.h>
#include "smoothstep.h"
#include <sys/time.h>
#include "cavityVolume.h"
#include "flowRate.h"
#include "massMatrix.h"
#include "setUpTraction.h"
#include "internalForce.h"
#include "enforceBC.h"
#include "mat2csv.h"
#include <math.h>
#include "internalForceBuckley.h"
#include "contactDetection.h"

#define REAL double

#define NUMTHREADS 4
/*  Function definitions */
jcv_point * readSites(int * numNodes, char fileName[64]);
int findInt(int * A, int val, int sizeArray);

int main(int argc, char** argv) {

	/*////////////////////////////////////////////////////////// */
	/*                                                           */
	/*			        Initialisation                           */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */


	struct timeval beginPre, endPre, endExp;

	gettimeofday(&beginPre, NULL);

	/* Initialsie random parameters  */
	FILE * fp;
	int i;
	char code;

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


	/*  Stretch rod */
	double dRodStop = 135; /*  130 mm  */
	double tStop = 0.5;
	double tRampRod = 1e-9;
	double vRod = 0;

	/*Delay pressure onset */
	double t_pressure = 1e-9;


	/*  Meshfree Paramters */
	double dMax = 2;
	double dMax_x = 2.5;
	double dMax_y = 2.5;	


	/*  Damping */
	double alpha = 0; 

	// stretch rod polynnomail
	double a0 = -2.2264e7;
	double a1 = 2.3704e7;
	double a2 = -9.3769e6;
	double a3 = 1.6212e6;
	double a4 =-9.7380e4;
	double a5 = -1.8801e3;
	double a6 = 559.3131;
	double a7 = 0.2565;



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


	/*  Stretch rod  */
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

	char opt[20] = "pDq0a2";
	char fileName[30] = "preform";
	double * points_out ;
	int * boundaryNodes;
	int numBoundary;
	int * nodalMarkers;
	int  numnodes;
	double * temperatures;
	trigen(&points_out,&boundaryNodes,opt,fileName,&numnodes,&numBoundary,&nodalMarkers,&temperatures);	

	fp = fopen("boundary.txt","w");
	for (int i = 0; i < numBoundary; ++i)
	{
		/* code */
		fprintf(fp,"%i\n",boundaryNodes[i]+1);
	}
	fclose(fp);

	// generate voronoi diagram
 	gpc_polygon ** voronoi = NULL;
 	callVoronoi(&voronoi,points_out,numnodes);
 	clipVoronoi(&voronoi,points_out,boundaryNodes,numnodes,numBoundary);



	/*  Set up efgBlock 
	 *  nodes
	 *  meshfree domain
	 *  */
	efg_block * efgBlock = malloc(1*sizeof(efg_block));
	strcpy(efgBlock->type,"radial");
	efgBlock->nodes = m_get(numnodes,2);
	for ( int i = 0 ; i < efgBlock->nodes->m ;i++){
		efgBlock->nodes->me[i][0] = points_out[2*i];
		efgBlock->nodes->me[i][1] = points_out[2*i+1];
	}
	/*  Set meshfree domain */
	efgBlock->di = v_get(numnodes);
	efgBlock->diX = v_get(numnodes);
	efgBlock->diY = v_get(numnodes);
	efgBlock->numnode = numnodes;
	printf("numnode = %i \n",efgBlock->numnode);
	efgBlock->dMax = v_get(3);
	efgBlock->dMax->ve[0] = dMax;
	efgBlock->dMax->ve[1] = dMax;
	efgBlock->dMax->ve[2] = dMax;
	setDomain(efgBlock);
	/*  Get traction and essential boundaries 
	 *  1 = traction
	 *  4 = essential
	 *  */


	/*  Traction boundary */
	boundaryLoad * traction = malloc(1*sizeof(boundaryLoad));
	traction->integrationType = 't';
	traction->type = 'p';
	traction->pressure = 0;
	traction->formulation ='e';
	getBoundary(&traction->nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,2);
	iv_foutput(stdout,traction->nodes);
	setUpTraction(traction,efgBlock);


	m_foutput(stdout,traction->coords);



	/*  Set up scni cells  */
	printf("#Set up SCNI cells\n");
	double site[2];
	SCNI * scni = malloc(efgBlock->numnode * sizeof(SCNI));
	for ( int i = 0 ; i < efgBlock->numnode ; i++){
		site[0] = efgBlock->nodes->me[i][0];
		site[1] = efgBlock->nodes->me[i][1];
		generateSCNI(voronoi[i],&scni[i],site,efgBlock);
	}



	/*  DIVERGENCE FREE CONDITION
	 *  For completlely interior nodes the sum of the B matrix should be 0 ( something like that )
	 *  */
	double totalArea = 0;
	int nodeNum = 33;


	printf("checking diveregence free condition at node %d, with corodinates %lf %lf \n", nodeNum, efgBlock->nodes->me[nodeNum][0],
		efgBlock->nodes->me[nodeNum][1]);
	MAT * divFree = m_get(4,2);
	for ( int i = 0 ; i < efgBlock->numnode ; i++){
		int res = findInt(scni[i].sfIndex->ive,nodeNum,scni[i].sfIndex->max_dim);
		if ( res != -1){
			divFree->me[0][0] += scni[i].area*scni[i].B->me[0][2*res];
			divFree->me[1][1] += scni[i].area*scni[i].B->me[1][2*res+1];
			divFree->me[2][0] += scni[i].area*scni[i].B->me[2][2*res];
			divFree->me[3][1] += scni[i].area*scni[i].B->me[3][2*res+1];
		}
		totalArea+= scni[i].area;	
	}

	printf("total area = %lf\n",totalArea);
	printf("#Divergence free condition\n");
	m_foutput(stdout,divFree);

	/*  Mass matrix */
	MAT * mass = m_get(efgBlock->numnode*2,efgBlock->numnode*2);
	massMatrix(mass,NULL,scni,rho, 'n', 'L', efgBlock->numnode);

	MAT * invMass = m_get(efgBlock->numnode*2,efgBlock->numnode*2);

	for ( int i = 0 ; i < invMass->m ; i++){
		invMass->me[i][i] = 1.0/mass->me[i][i];
	}







	/*  Set up essential boundary  */
	EBC * eb1 = malloc(1*sizeof(EBC));
	eb1->dofFixed = 3;
	getBoundary(&eb1->nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,5);
	iv_addNode(eb1->nodes,traction->nodes->ive[0 ],'s');
	setUpBC(eb1,mass,efgBlock);
	m_foutput(stdout,eb1->coords);

	/*  Set up essential boundary  */
	EBC * eb2 = malloc(1*sizeof(EBC));
	eb2->dofFixed = 1;
	getBoundary(&eb2->nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,4);
	iv_addNode(eb2->nodes,traction->nodes->ive[traction->nodes->max_dim - 1],'e');
	setUpBC(eb2,mass,efgBlock);

	m_foutput(stdout,eb2->coords);

	/*  Set up essential boundary  */
	int numB3 = 15;
	EBC * eb3 = malloc(1*sizeof(EBC));
	eb3->dofFixed = 2;
	eb3->nodes = iv_get(numB3);
	for ( int i = 0 ; i < eb3->nodes->max_dim ; i++){
		eb3->nodes->ive[i] = traction->nodes->ive[traction->nodes->max_dim -1 - i];
	}
	setUpBC(eb3,mass,efgBlock);

	m_foutput(stdout,eb3->coords);
	MAT * testPoint = m_get(1,2);



	MAT * point = m_get(1,2);
	EFG_SF * phiContact = malloc(eb3->nodes->max_dim * sizeof(EFG_SF));

	define_support(phiContact,eb3->coords,efgBlock);
	MLS_shapefunction(eb3->coords,phiContact,efgBlock);
	/*  Get transformation matrix */
	/* Lambda  */

	MAT * Lambda = m_get(2*efgBlock->numnode,2*efgBlock->numnode);
	EFG_SF * _phi = malloc(1*sizeof(EFG_SF));
	for ( int i = 0 ; i < efgBlock->numnode ; i++){
		point->me[0][0] = efgBlock->nodes->me[i][0];
		point->me[0][1] = efgBlock->nodes->me[i][1];
		define_support(_phi,point,efgBlock);
		MLS_shapefunction(point,_phi,efgBlock);
		for ( int k = 0; k < _phi->index->max_dim ; k++){
			Lambda->me[2*i][2*_phi->index->ive[k]] += _phi->phi->ve[k];
			Lambda->me[2*i+1][2*_phi->index->ive[k]+1] += _phi->phi->ve[k];
		}
	}	


	/*  Timing */
	gettimeofday(&endPre, NULL);
	long elapsedPre = (endPre.tv_sec-beginPre.tv_sec) + (endPre.tv_usec-beginPre.tv_usec)/1000000.0;

	///////////////////////////////////////////////////////////////

	/*////////////////////////////////////////////////////////// */
	/*                                                           */
	/*			      EXPLICIT TIME STEPPNG                      */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */

	/*  Initialise variables */
	timestep * timeStep = malloc(1*sizeof(timestep));
	timeStep->deltaT = 1e-6;
	timeStep->tmax = 1;
	/*  time step  */
	double t_n_1;
	double t_n = 0;
	double t_n_h;

	/*  Kinematic variables */
	VEC * a_n_1 = v_get(efgBlock->numnode*2);
	VEC * a_n = v_get(efgBlock->numnode*2);
	VEC * d_n_1 = v_get(efgBlock->numnode*2);
	VEC * d_n = v_get(efgBlock->numnode*2);
	VEC * v_n_1 = v_get(efgBlock->numnode*2);
	VEC * v_n_h = v_get(efgBlock->numnode*2);

	VEC * v_n_mh = v_get(efgBlock->numnode*2);
	VEC * v_n = v_get(efgBlock->numnode*2);
	VEC * deltaDisp = v_get(efgBlock->numnode*2);
	VEC * nodalDisp = v_get(efgBlock->numnode*2);
	MAT * updatedNodes = m_get(efgBlock->numnode,2);

	MAT * updatedNodesT = m_get(efgBlock->numnode,2);

	/*  Force variables */
	VEC * Fext_n_1 = v_get(efgBlock->numnode*2);
	VEC * Fext_n = v_get(efgBlock->numnode*2);
	VEC * Fint_n_1 = v_get(efgBlock->numnode*2);
	VEC * Fint_n = v_get(efgBlock->numnode*2);
	VEC * Fnet = v_get(efgBlock->numnode*2);
	VEC * Fdamp = v_get(efgBlock->numnode*2);

	/* Boundary conditions */
	eb1->uBar1 = v_get(eb1->nodes->max_dim);
	eb1->uBar2 = v_get(eb1->nodes->max_dim);
	eb2->uBar1 = v_get(eb2->nodes->max_dim);
	eb2->uBar1 = v_get(eb2->nodes->max_dim);
	eb3->uBar2 = v_get(eb3->nodes->max_dim);
	eb3->uBar2 = v_get(eb3->nodes->max_dim);

	VEC * v_correct = v_get(efgBlock->numnode*2);

	/*  Energy  */
	double Wkin_n_1;
	double Wext_n_1;
	double Wint_n_1;
	double Wint_n;
	double Wext_n;
	double Wbal;

	/*  Iteration counter */
	int n= 0;
	/*  File write counter */
	int writeFreq = 100;
	int fileCounter = 1;



	/*  Cavity pressure paramters */
	double chokedMassRate =  flowRate(0.528, tLine, pLine*1e6, rLine, aReduced, gammaLine);
	double volumeInitial = cavityVolume(traction->nodes,efgBlock->nodes);
	double massInitial = (P0*volumeInitial)/(rLine*1000*tLine);
	double massAir = massInitial;
	double pRatio;
	double pre_n = 0;
	double pre_n_1 = 0; 
	double volume = 0;
	double volume_t = 0;

	/*  For writing to file */
	fp = fopen("pressureTime.txt","w");
	fprintf(fp,"%lf %lf\n",0.0,0.0);
	fclose(fp);


	/*  Contact with stretch rod */
	double dispRod = 0;
	double dispRodStop = 0;
	double distanceProj = 0;
	MAT * msNormal = m_get(1,2);
	double f1Cor = 0;
	double f2Cor = 0;

	double beta_1 = 10;
	double beta_2 = 10;


	/*  Buckley state  */

	State * stateOld = malloc(efgBlock->numnode *sizeof *stateOld  );


	for ( int i = 0 ; i < efgBlock->numnode; i++){
		stateOld[i] =  malloc(sizeof *stateOld[i] );
		stateOld[i]->Fbar = m_get(3,3);
		m_ident(stateOld[i]->Fbar); 
		stateOld[i]->Dbar = m_get(3,3);
		stateOld[i]->Bbar = m_get(3,3);
		m_ident(stateOld[i]->Bbar);
		stateOld[i]->Wbar = m_get(3,3);
		stateOld[i]->Sc = m_get(3,3);
		stateOld[i]->Sb = m_get(3,3);
		stateOld[i]->eigValDBar = v_get(3);
		stateOld[i]->eigValDBar_p = v_get(3);
		stateOld[i]->trace_D = 0;
		stateOld[i]->critLambdaBar = 100.00;
		stateOld[i]->mSigma = 0;
		stateOld[i]->lambdaBar = v_get(3);
		stateOld[i]->lambdaNMax = 1;
		stateOld[i]->temperature = temperatures[i] ;
		stateOld[i]->ep_true = m_get(3,3);
		printf("temperatures  of node %d = %lf\n",i,temperatures[i]);
		stateOld[i]->Jacobian = 1;
		stateOld[i]->F = m_get(3,3);
		stateOld[i]->eigVecDBar = m_get(3,3);
		m_ident(stateOld[i]->eigVecDBar);
		m_ident(stateOld[i]->F);
		stateOld[i]->D = m_get(3,3);
		stateOld[i]->W = m_get(3,3);


		stateOld[i]->Drot = m_get(3,3);
		stateOld[i]->R = m_get(3,3);
		m_ident(stateOld[i]->R);


	}

	mv_mlt(Lambda,d_n_1,nodalDisp);
	for ( int i = 0 ; i < efgBlock->numnode ; i++){
		updatedNodes->me[i][0]  = efgBlock->nodes->me[i][0] + nodalDisp->ve[2*i] ;
		updatedNodes->me[i][1]  = efgBlock->nodes->me[i][1] + nodalDisp->ve[2*i+1] ;	

		updatedNodesT->me[i][0]  = efgBlock->nodes->me[i][0] + nodalDisp->ve[2*i] ;
		updatedNodesT->me[i][1]  = efgBlock->nodes->me[i][1] + nodalDisp->ve[2*i+1] ;	

	}

	double disp_rod = 0;
	/*  Explicit Loop */
	while ( t_n < tMax){

		/*  Update time step */
		t_n_1 = t_n + deltaT;
		t_n_h = (1.00/2.00)*(t_n + t_n_1);

		/*  Make a time step  */ 
		v_mltadd(v_n,a_n,(double)0.5*deltaT,v_n_h);
		v_mltadd(d_n,v_n_h,deltaT,d_n_1);

		if ( disp_rod < 132){
		/*  Update stretch rod */
			disp_rod = a0*pow(t_n_1,7) + a1*pow(t_n_1,6) + a2*pow(t_n_1,5) + a3*pow(t_n_1,4) + a4*pow(t_n_1,3) + a5*pow(t_n_1,2) +
			a6*pow(t_n_1,1) + a7;
			//disp_rod = vRod*t_n_1;

			for ( int i = 0 ; i < srNodes->m ; i++){
	
				srNodes->me[i][1] = srNodes_O->me[i][1] - disp_rod;

			}
			dispRod = disp_rod;
		}
		/*  Stretch rod boundary conditions */

		for ( int i = 0 ; i < eb3->nodes->max_dim ; i++){


			testPoint->me[0][0] = updatedNodes->me[eb3->nodes->ive[i]][0];
			testPoint->me[0][1] = updatedNodes->me[eb3->nodes->ive[i]][1];

			distanceProj = contactDetection(testPoint,srNodes,msNormal);

			//printf("distance = %lf with normal {%lf,%lf}\n",distanceProj,msNormal->me[0][0],msNormal->me[0][1]);


			if (distanceProj > 0){

				f1Cor = 1*(2*distanceProj*msNormal->me[0][0]*mass->me[eb3->nodes->ive[i]*2][eb3->nodes->ive[i]*2])/pow(deltaT,2);
				f2Cor = 1*(2*distanceProj*msNormal->me[0][1]*mass->me[eb3->nodes->ive[i]*2][eb3->nodes->ive[i]*2])/pow(deltaT,2);



				for ( int k = 0 ; k < phiContact[i].index->max_dim ; k++){
					Fnet->ve[2*phiContact[i].index->ive[k]] += phiContact[i].phi->ve[k]*f1Cor; 
					Fnet->ve[2*phiContact[i].index->ive[k]+1] += phiContact[i].phi->ve[k]*f2Cor; 
				}

			}


				// penalty contact
				//double gamma_N = vRod - v_n_h->ve[eb3->nodes->ive[i]*2];
				//double g_N = distanceProj;
				// penalty method
				//double p_bar = g_N gamma_N*beta_2*heaviside(gamma_N)


		}
		mv_mlt(invMass,Fnet,a_n);
		/*  Find a corrective acceleration - method in pronto 3D manual*/
		/*  Make a time step  */ 
		v_mltadd(v_n,a_n,(double)0.5*deltaT,v_n_h);
		v_mltadd(d_n,v_n_h,deltaT,d_n_1);


		/*  Implement BCs */
		enforceBC(eb1,d_n_1); 
		// find velocity correction
		sv_mlt(1.00/(2*deltaT),eb1->uCorrect1,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k] += v_correct->ve[k];
		}

		sv_mlt(1.000/(2*deltaT),eb1->uCorrect2,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k+1] += v_correct->ve[k];
		}


		// Symmetry boundary /
		enforceBC(eb2,d_n_1); 
		sv_mlt(1.00/(2*deltaT),eb2->uCorrect1,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k] += v_correct->ve[k];
		}


		/*  Update nodal poistions */
		mv_mlt(Lambda,d_n_1,nodalDisp);
		for ( int i = 0 ; i < efgBlock->numnode ; i++){
			updatedNodes->me[i][0]  = efgBlock->nodes->me[i][0] + nodalDisp->ve[2*i] ;
			updatedNodes->me[i][1]  = efgBlock->nodes->me[i][1] + nodalDisp->ve[2*i+1] ;

			if ( n % 10 == 0 )
				{
				updatedNodesT->me[i][0]  = efgBlock->nodes->me[i][0] + nodalDisp->ve[2*i] ;
				updatedNodesT->me[i][1]  = efgBlock->nodes->me[i][1] + nodalDisp->ve[2*i+1] ;	
			}

		}

		// if ( n % 100000 == 0) 
		// {	
		// 	 gpc_polygon ** voronoi_1;
		// 	 callVoronoi(&voronoi_1,updatedNodes->base,numnodes);
 	// 		 clipVoronoi(&voronoi_1,updatedNodes->base,boundaryNodes,numnodes,numBoundary);
		// }


		/*  Find Cavity volume */
		volume = cavityVolume(traction->nodes,updatedNodes);

		if ( t_n_1 < t_pressure)
		{
			pre_n_1 = 0;
		}else{


		/*  Find Cavity pressure */
		pRatio = pre_n/pLine;
		if ( pRatio <= 0.528){
			massAir += chokedMassRate*deltaT;
		}else{
			massAir += flowRate(pRatio,tLine,pLine*1e6, rLine, aReduced,gammaLine)*deltaT;
		}
		pre_n_1 = ((P0*(volume - volumeInitial) + 1000*massAir*rLine*tLine)/(volume+vDead));
		}


		//pre_n_1 = 0.5*smoothstep(t_n_1,0.2,t_pressure);

		/*  Update pressure load */
		traction->pressure = -pre_n_1;
		externalForce(Fext_n_1,traction,scni,updatedNodes);


		/*  Internal force */
		//internalForce(Fint_n_1,&scni,d_n_1,matParams,efgBlock->numnode);
		internalForceBuckley(Fint_n_1,scni,d_n_1,matParams,critLambdaParams,stateOld,deltaT,efgBlock->numnode,t_n_1);

		/*  Damping */
		mv_mlt(mass,v_n_h,Fdamp);
		sv_mlt(alpha,Fdamp,Fdamp);


		/*  Balance of forces */
		v_sub(Fext_n_1,Fint_n_1,Fnet);
		v_sub(Fnet,Fdamp,Fnet);

		/*  Find acceleration */
		mv_mlt(invMass,Fnet,a_n_1);

		/*  Integer time step velocity */
		v_mltadd(v_n_h,a_n_1,t_n_1-t_n_h,v_n_1);	


		// update nodal positions
		if ( n % writeFreq == 0 ){
			char * filename[50];
			snprintf(filename, 50, "displacement_%d%s",fileCounter,".txt");
			mat2csv(updatedNodes,"./Displacement",filename);

			snprintf(filename, 50, "srRod_%d%s",fileCounter,".txt");
			mat2csv(srNodes,"./srRod",filename);

			fp = fopen("pressureTime.txt","a");
			fprintf(fp,"%lf %lf\n",t_n_1,pre_n_1);
			fclose(fp);

			fileCounter++;
		}


		/*  Find energy */
		v_sub(d_n_1,d_n,deltaDisp);
		// External Energy
		Wext_n_1 = Wext_n + 0.5*(in_prod(deltaDisp,Fext_n) + in_prod(deltaDisp,Fext_n_1));
		// Internl Energy
		Wint_n_1 = Wint_n + 0.5*(in_prod(deltaDisp,Fint_n) + in_prod(deltaDisp,Fint_n_1));
		// Kinetic Energy = 1/2 v^T * M * v
		Wkin_n_1 = 0 ;
		for (int i = 0 ; i < 2* efgBlock->numnode ; i++){

			Wkin_n_1 = Wkin_n_1 + 0.5*(v_n_1->ve[i]*mass->me[i][i]*v_n_1->ve[i]);
		}
		// Find the energy balance.
		Wbal = fabs(Wkin_n_1 + Wint_n_1 - Wext_n_1) ;
		// Find maximum of the energies, used to normalise the balance. Uses the max macro defined in matrix.h
		double W_max = max(Wext_n_1,Wint_n_1);
		W_max = max(W_max,Wkin_n_1);
		/*  Update counters */
		t_n = t_n_1;
		// Store energy accumulated
		Wext_n = Wext_n_1;
		Wint_n = Wint_n_1;
		/*  Store previous volume */
		volume_t = volume; 
		// Store previous time step quanities for the kinematic, and force variables.
		v_copy(Fint_n_1,Fint_n);
		v_copy(Fext_n_1,Fext_n);
		v_copy(v_n_h,v_n_mh);
		v_copy(d_n_1,d_n);
		v_copy(v_n_1,v_n);
		v_copy(a_n_1,a_n);
		pre_n = pre_n_1;
		// update iteration counter
		n++	;
		printf("%i  \t  %lf %10.2E %lf %lf %lf \n",n,t_n,Wbal/W_max,pre_n_1,dispRod,volume/1e3);




	}


	/*////////////////////////////////////////////////////////// */
	/*////////////////////////////////////////////////////////// */


	gettimeofday(&endExp, NULL);

	long elapsedExp = (endExp.tv_sec-endPre.tv_sec) + (endExp.tv_usec-endPre.tv_usec)/1000000.0;


	printf("Post processing took %ld seconds \n Explicit routine took %ld seconds\n",elapsedPre,elapsedExp);




	// /*////////////////////////////////////////////////////////// */
	// /*			      POSTPROCESSING STAGE                       */
	// /*                                                           */
	// /*////////////////////////////////////////////////////////// */







	/*////////////////////////////////////////////////////////// */
	/*////////////////////////////////////////////////////////// */



	/*  Program exit */
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