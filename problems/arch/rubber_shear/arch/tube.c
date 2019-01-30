#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define JC_VORONOI_IMPLEMENTATION
#include "structure.h"
#include "vGeom.h"
#include "triangle.h"
#include "setDomain.h"
#include "tri_mesh.h"
#include "matrix.h"
#include "matrix2.h"
#include "generateSCNI.h"
#include "getBoundary.h"
#include "iv_addNode.h"
#include "setUpBC.h"
#include "jc_voronoi.h"
#include "clipVoronoi.h"
#include "externalForce.h"
#include <pthread.h>
#include "smoothstep.h"
#include <sys/time.h>
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
	jcv_diagram diagram;
	jcv_graphedge* graph_edge;
	memset(&diagram, 0, sizeof(jcv_diagram));
	char code;
	jcv_point  pInt ; 
	int numBoundary;
	int * boundaryNodes ;

	/*  Material parameters */
	const double rho = 1000e-9;
		/*  Material struct */
	VEC * matParams = v_get(4);
	matParams->ve[0] = 0.373;
	matParams->ve[1] = -0.031;
	matParams->ve[2] = 0.005;
	matParams->ve[3] = 1e5;

	/*  Time step variables */

	double deltaT = 1e-6;
	double tMax = 10;

	double tRamp =3;
	double tRamp2 = 1.5;


	/*  Meshfree Paramters */
	double dMax = 3;
	double dMax_x = 2.5;
	double dMax_y = 2.5;	


	/*  Loading parameters */
	double pressure = -0.2;
	double pre1 = -0.2;
	double pre2 = -0.3;
	/*  Plotting parameters */


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

	struct triangulateio * out = malloc(1*sizeof(struct triangulateio )); ;

	/*  Manual triangulation */
	double rin = 60;
	double rout = 80;
	double width = 2;
	int numPoint = 15;

	out->pointlist = malloc(2*2*numPoint*sizeof(double));
	out->pointmarkerlist = malloc(2*numPoint*sizeof(int));
	out->segmentlist = malloc(2*2*numPoint*sizeof(int));
	boundaryNodes = malloc(2*numPoint*sizeof(int));
	out->numberofpoints = 2*numPoint ;
	numBoundary = 2*numPoint;
	for ( i = 0 ; i < numPoint*2 ; i++){
		if ( i < numPoint ){

		out->pointlist[2*i] = rin + (i)*(rout-rin)/(numPoint-1);
		out->pointlist[2*i+1] = width;
		out->pointmarkerlist[i] = 3;
		}else{
		out->pointlist[2*i] = rout - (i-numPoint)*(rout-rin)/(numPoint-1);
		out->pointlist[2*i+1] =0 ;
		out->pointmarkerlist[i] = 4;

		} 
		boundaryNodes[i] = i;
	}
	out->pointmarkerlist[0] = 2;
	out->pointmarkerlist[2*numPoint - 1] = 2;



	jcv_point * points = malloc(out->numberofpoints * sizeof(jcv_point));

	for ( int i = 0 ; i < out->numberofpoints ; i++){
		points[i].x = out->pointlist[2*i];
		points[i].y = out->pointlist[2*i+1];
	}
	int numNodes = out->numberofpoints;


	for ( int i = 0 ; i < numNodes ; i++){
		printf("%lf\t%lf\n",points[i].x,points[i].y);
	}

	FILE * ptr;
	ptr = fopen("boundary.txt","w");
	for (int i = 0 ; i < numBoundary ; i++){
		fprintf(ptr,"%i\n",boundaryNodes[i] + 1);

	}
	fclose(ptr);



	/*  Set up efgBlock 
	 *  nodes
	 *  meshfree domain
	 *  */
	efg_block * efgBlock = malloc(1*sizeof(efg_block));
	strcpy(efgBlock->type,"radial");
	efgBlock->nodes = m_get(numNodes,2);
	for ( int i = 0 ; i < efgBlock->nodes->m ;i++){
		efgBlock->nodes->me[i][0] = points[i].x;
		efgBlock->nodes->me[i][1] = points[i].y;
	}
	/*  Set meshfree domain */
	efgBlock->di = v_get(numNodes);
	efgBlock->diX = v_get(numNodes);
	efgBlock->diY = v_get(numNodes);
	efgBlock->numnode = numNodes;
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
	traction->nodes = iv_get(2);
	traction->nodes->ive[0] = 2*numPoint - 1 ;
	traction->nodes->ive[1] = 0;
	iv_foutput(stdout,traction->nodes);
	setUpTraction(traction,efgBlock);
	m_foutput(stdout,traction->coords);
	/*  Generate Voronoi diagram */
	printf("genereating voronoi diagram\n");
	jcv_diagram_generate(numNodes, points,(jcv_rect *) NULL, &diagram);


	/*  Get sites */
	jcv_site* sites = malloc(numNodes*sizeof(jcv_site)) ;
	sites = jcv_diagram_get_sites(&diagram);

	/*  Clipping voronoi */
	printf("clipping voronoi \n");
	/*  Clip with boundaryNodes */
	clipVoronoi(sites,points,boundaryNodes,numBoundary,numNodes);


	/*  Write edges */


	fp = fopen("simple.edges","w");
	printf("# Write Edges to File .edges \n");
	sites = jcv_diagram_get_sites(&diagram);
	for (i=0; i<diagram.numsites; i++) {

		graph_edge = sites[i].edges;
		while (graph_edge) {
			// This approach will potentially print shared edges twice
			fprintf(fp,"%f %f\n", (double)graph_edge->pos[0].x, (double)graph_edge->pos[0].y);
			fprintf(fp,"%f %f\n", (double)graph_edge->pos[1].x, (double)graph_edge->pos[1].y);
			graph_edge = graph_edge->next;
		}
		fprintf(fp,"%s\n","\\");
	}

	fclose(fp);


	// write sites to file
	fp = fopen("simple.sites","w");
	printf("# Write Seeds to File .sites \n");
	for (i=0; i<numNodes; i++) {
		fprintf(fp,"%f,%f\n", (double)points[i].x, (double)points[i].y);
	}
	/*  Set up traction boundary */

	fclose(fp);

	/*  Set up scni cells  */
	printf("#Set up SCNI cells\n");
	SCNI * scni = malloc(efgBlock->numnode * sizeof(SCNI));
	for ( int i = 0 ; i < efgBlock->numnode ; i++){
		generateSCNI(&sites[i],&scni[i],efgBlock);
	}



	/*  DIVERGENCE FREE CONDITION
	 *  For completlely interior nodes the sum of the B matrix should be 0 ( something like that )
	 *  */
	double totalArea = 0;
	int nodeNum = 121;
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
	eb1->dofFixed = 1;
	eb1->nodes = iv_get(2);
	eb1->nodes->ive[0] = 2*numPoint - 1 ;
	eb1->nodes->ive[1] = 0;
	setUpBC(eb1,mass,efgBlock);
	m_foutput(stdout,eb1->coords);



	/*  Set up essential boundary  */
	EBC * eb2 = malloc(1*sizeof(EBC));
	eb2->dofFixed = 2;
	getBoundary(&eb2->nodes,boundaryNodes,numBoundary,out,3);
	iv_addNode(eb2->nodes,eb1->nodes->ive[eb1->nodes->max_dim -1 ],'s');
	setUpBC(eb2,mass,efgBlock);
	m_foutput(stdout,eb2->coords);

	/*  Set up essential boundary  */
	EBC * eb3 = malloc(1*sizeof(EBC));
	eb3->dofFixed = 2;
	getBoundary(&eb3->nodes,boundaryNodes,numBoundary,out,4);
	iv_addNode(eb3->nodes,eb1->nodes->ive[0],'e');
	setUpBC(eb3,mass,efgBlock);

	m_foutput(stdout,eb3->coords);

	int nodeSearch = 0;
	for ( int i = 0 ; i < efgBlock->numnode ; i++){
		double tempX = scni[i].coords[0];

		if ( tempX == rin){
			nodeSearch = i;
		}
	}
	printf("nodeSearch = %i\n",nodeSearch);

	/*  Get transformation matrix */
	/* Lambda  */

	MAT * Lambda = m_get(2*efgBlock->numnode,2*efgBlock->numnode);
	EFG_SF * _phi = malloc(1*sizeof(EFG_SF));
	MAT * point = m_get(1,2);
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


	/*  Force variables */
	VEC * Fext_n_1 = v_get(efgBlock->numnode*2);
	VEC * Fext_n = v_get(efgBlock->numnode*2);
	VEC * Fint_n_1 = v_get(efgBlock->numnode*2);
	VEC * Fint_n = v_get(efgBlock->numnode*2);
	VEC * Fnet = v_get(efgBlock->numnode*2);

	VEC * stressVoigt = v_get(5);
	/* Boundary conditions */
	eb1->uBar1 = v_get(eb1->nodes->max_dim);
	eb1->uBar2 = v_get(eb1->nodes->max_dim);
	eb2->uBar1 = v_get(eb2->nodes->max_dim);
	eb2->uBar2 = v_get(eb2->nodes->max_dim);
	eb3->uBar1 = v_get(eb3->nodes->max_dim);
	eb3->uBar2 = v_get(eb3->nodes->max_dim);

	VEC * v_correct = v_get(efgBlock->numnode*2);

	/*  Energy  */
	double Wkin_n_1;
	double Wext_n_1;
	double Wint_n_1;
	double Wint_n;
	double Wext_n;
	double Wbal;

	double tStop = 0;
	double preStop = 0;
	/*  Iteration counter */
	int n= 0;
	/*  File write counter */
	int writeFreq = 100;
	int fileCounter = 1;

	double pre_n_1 = 0; 

	
	/*  For writing to file */
	fp = fopen("loadDisp.txt","w");
	fprintf(fp,"%lf %lf\n",0,0);
	fclose(fp);

	while ( t_n < tMax){

		/*  Update time step */
		t_n_1 = t_n + deltaT;
		t_n_h = (1.00/2.00)*(t_n + t_n_1);

		/*  Make a time step  */ 
		v_mltadd(v_n,a_n,(double)0.5*deltaT,v_n_h);
		v_mltadd(d_n,v_n_h,deltaT,d_n_1);


	
		/*  Implement BCs */
		enforceBC(eb2,d_n_1); 
		// find velocity correction
		sv_mlt(1.00/(2*deltaT),eb2->uCorrect2,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k+1] += v_correct->ve[k];
		}

		enforceBC(eb3,d_n_1); 
		sv_mlt(1.000/(2*deltaT),eb3->uCorrect2,v_correct);
		for ( int k = 0 ; k < v_correct->max_dim; k++){
			v_n_h->ve[2*k+1] += v_correct->ve[k];
		}

		/*  Update nodal positions */
		
		mv_mlt(Lambda,d_n_1,nodalDisp);
		for ( int i = 0 ; i < efgBlock->numnode ; i++){
			updatedNodes->me[i][0]  = efgBlock->nodes->me[i][0] + nodalDisp->ve[2*i] ;
			updatedNodes->me[i][1]  = efgBlock->nodes->me[i][1] + nodalDisp->ve[2*i+1] ;
		}





		/*  Update tip load */

		if ( nodalDisp->ve[2*eb1->nodes->ive[0]] < 80){
			pre_n_1 = pre1 * smoothstep(t_n_1,tRamp,0);
			tStop = t_n_1;
			preStop = pre_n_1;
		}
		else{
			pre_n_1 = pre2 * smoothstep(t_n_1,0.5+tStop,tStop+deltaT) + preStop;
		}
		traction->pressure = pre_n_1;
		externalForce(Fext_n_1,traction,scni,updatedNodes);

		
		internalForce(Fint_n_1,scni,d_n_1,matParams,efgBlock->numnode);
	
		
	
		

		/*  Balance of forces */
		v_sub(Fext_n_1,Fint_n_1,Fnet);
		
		/*  Find acceleration */
		mv_mlt(invMass,Fnet,a_n_1);

		/*  Integer time step velocity */
		v_mltadd(v_n_h,a_n_1,t_n_1-t_n_h,v_n_1);	


		// update nodal positions
		if ( n % writeFreq == 0 ){
			char * filename[50];
			snprintf(filename, 50, "displacement_%d%s",fileCounter,".txt");
			mat2csv(updatedNodes,"./Displacement",filename);
			
			fp = fopen("loadDisp.txt","a");
			fprintf(fp,"%lf %lf\n",nodalDisp->ve[2*eb1->nodes->ive[0]],-pre_n_1);
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
		// Store previous time step quanities for the kinematic, and force variables.
		v_copy(Fint_n_1,Fint_n);
		v_copy(Fext_n_1,Fext_n);
		v_copy(v_n_h,v_n_mh);
		v_copy(d_n_1,d_n);
		v_copy(v_n_1,v_n);
		v_copy(a_n_1,a_n);
		// update iteration counter
		n++	;
		printf("%i  \t  %lf %lf \t  %10.2E \n",n,t_n,pre_n_1,Wbal/W_max);

	}


	/*////////////////////////////////////////////////////////// */
	/*////////////////////////////////////////////////////////// */


	gettimeofday(&endExp, NULL);

	long elapsedExp = (endExp.tv_sec-endPre.tv_sec) + (endExp.tv_usec-endPre.tv_usec)/1000000.0;


	printf("Post processing took %ld seconds \n Explicit routine took %ld seconds\n",elapsedPre,elapsedExp);




	/*////////////////////////////////////////////////////////// */
	/*			      POSTPROCESSING STAGE                       */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */

	
	/*  Free allocated memory */
	jcv_diagram_free(&diagram);






	/*////////////////////////////////////////////////////////// */
	/*////////////////////////////////////////////////////////// */



	/*  Program exit */
	return 0;
}
