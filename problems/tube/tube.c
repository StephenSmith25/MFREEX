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
#include "Force/Internal/internalForce_Buckley.h"
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
#include "Material/Buckley/new_Buckley_State.h"
#include "Boundary/Contact/contactDetection.h"
#include "Boundary/Displacement/setUpBC.h"
#include "Boundary/Displacement/enforceBC.h"
#include "Integration/SCNI/scni_update_B.h"


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
	int * nodalMarkers = out->pointmarkerlist;

	int numnodes = out->numberofpoints;
	int dim = 2;
	int is_AXI = 1;



	MAT * xI = m_get(numnodes,dim);
	// create nodes
	for ( int i = 0 ; i < numnodes ; i++)
	{
		xI->me[i][0] = out->pointlist[2*i];
		xI->me[i][1] = out->pointlist[2*i+1];

	}

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
	double dmax = 2;
	int constant_support_size = 1;
	VEC * dI = v_get(xI->m);

	// meshfree domain
	meshfreeDomain mfree = {.nodes = xI, .di = dI, .num_nodes = xI->m, .dim = dim, .IS_AXI = is_AXI};
	setDomain(&mfree,constant_support_size, dmax);

	v_foutput(stdout,mfree.di);

	

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

		//nodal_mass->ve[i] = volume*rho;
		//inv_nodal_mass->ve[i] = 1.000/nodal_mass->ve[i];

	}
	for ( int i = 0 ; i < mfree.num_nodes ; i++)
	{
		inv_nodal_mass->ve[i] = 1.00/nodal_mass->ve[i];
	}
	printf("total mass = %lf (g) \n", v_sum(nodal_mass)*1000);
	v_foutput(stdout,nodal_mass);
	


	/*  Set up essential boundary  */
	EBC * eb1 = malloc(1*sizeof(EBC));
	eb1->dofFixed = 1;
	eb1->nodes = iv_get(2);
	eb1->nodes->ive[0] = 2*numPoint - 1 ;
	eb1->nodes->ive[1] = 0;
	int num_nodes_eb1 = eb1->nodes->max_dim;
	setUpBC(eb1,inv_nodal_mass,&mfree);
	m_foutput(stdout,eb1->coords);



	/*  Set up essential boundary  */
	EBC * eb2 = malloc(1*sizeof(EBC));
	eb2->dofFixed = 2;
	getBoundary(&eb2->nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,3);
	iv_addNode(eb2->nodes,eb1->nodes->ive[eb1->nodes->max_dim -1 ],'s');
	setUpBC(eb2,inv_nodal_mass,&mfree);
	m_foutput(stdout,eb2->coords);

	/*  Set up essential boundary  */
	EBC * eb3 = malloc(1*sizeof(EBC));
	eb3->dofFixed = 2;
	getBoundary(&eb3->nodes,boundaryNodes,numBoundary,nodalMarkers,numnodes,4);
	iv_addNode(eb3->nodes,eb1->nodes->ive[0],'e');
	setUpBC(eb3,inv_nodal_mass,&mfree);

	m_foutput(stdout,eb3->coords);




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

	free_shapefunction_container(sf_nodes);


	gettimeofday(&endPre, NULL);
	long elapsedPre = (endPre.tv_sec-beginPre.tv_sec) + (endPre.tv_usec-beginPre.tv_usec)/1000000.0;

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
	VEC * deltaDisp = v_get(num_dof);
	VEC * nodalDisp = v_get(num_dof);
	MAT * updatedNodes = m_get(num_dof,dim);


	/*  Force variables */
	VEC * Fext_n_1 = v_get(num_dof);
	VEC * Fext_n = v_get(num_dof);
	VEC * Fint_n_1 = v_get(num_dof);
	VEC * Fint_n = v_get(num_dof);
	VEC * Fnet = v_get(num_dof);

	VEC * stressVoigt = v_get(5);
	/* Boundary conditions */
	eb1->uBar1 = v_get(eb1->nodes->max_dim);
	eb1->uBar2 = v_get(eb1->nodes->max_dim);
	eb2->uBar1 = v_get(eb2->nodes->max_dim);
	eb2->uBar2 = v_get(eb2->nodes->max_dim);
	eb3->uBar1 = v_get(eb3->nodes->max_dim);
	eb3->uBar2 = v_get(eb3->nodes->max_dim);

	VEC * v_correct = v_get(num_dof);

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
	fprintf(fp,"%lf %lf\n",0.00,0.00);
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





		/*  Update pressure load */

		if ( nodalDisp->ve[2*eb1->nodes->ive[0]] < 80){
			pre_n_1 = pre1 * smoothstep(t_n_1,tRamp,0);
			tStop = t_n_1;
			preStop = pre_n_1;
		}
		else{
			pre_n_1 = pre2 * smoothstep(t_n_1,0.5+tStop,tStop+deltaT) + preStop;
		}


		
		internalForce(Fint_n_1,scni,d_n_1,matParams,efgBlock->numnode);
	
		
		update_pressure_boundary(pB, updatedNodes);
		assemble_pressure_load(Fext_n_1, pre_n_1, pB);
		

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



	/*////////////////////////////////////////////////////////// */
	/*////////////////////////////////////////////////////////// */



	/*  Program exit */
	return 0;
}
