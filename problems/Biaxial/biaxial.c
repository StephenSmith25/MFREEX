#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "Material/Buckley/buckleyStress.h"
#include "Material/material.h"
#include <sys/time.h>
#include "Material/Hyperelastic/yeoh.h"
#include "Deformation/update_Polar_Decomposition.h"
#include "Deformation/rotate_tensor.h"
#include "Material/Buckley/buckleyStress.h"
#include "matrix.h"
#include "matrix2.h"
#include <omp.h>
#include "Integration/SCNI/generate_scni.h"
#include "Material/material.h"
#include "mat2csv.h"
#include "m_inverse_small.h"
#include "Material/material.h"
#include "Integration/gMat.h"
#include "Deformation/velocity_grad.h"
#include "symmeig_small.h"
#include "dsyevh3.h"
#include "dsyevv3.h"
#include "dsyevq3.h"

// material
const char * MATERIAL = "BUCKLEY";
const int BUCKLEY_MATERIAL = 1;
const int PLASTIC_MATERIAL = 0;

// deformaton
const int SR = 4;
const double TEMPERATURE = 85;
char * DEFORMATION_MODE = "SIMPLE_SHEAR";
const int DIM = 3;
const int IS_AXI = 0;

// time step
const double DT = 1e-5;
const double TMAX = 10;
const double PEAK_STRAIN = 5;


int main(void)
{



	/*  Material struct */
	VEC * matParams = v_get(31);
	matParams->ve[0] = 2.814e-3; // VS
	matParams->ve[1] = 0.526e-3; // VP
	matParams->ve[2] = 1.71e6; // mu*_0
	matParams->ve[3] = 328.76; // Tinf
	matParams->ve[4] = 358.15; // T*
	matParams->ve[5] = matParams->ve[4]; // Tf*
	matParams->ve[6] = (67.47); // Cv
	matParams->ve[7] = 1.23e5; // H0
	matParams->ve[8] = 8.314; // R
	matParams->ve[9] = 1.8e9; // Kb
	matParams->ve[10] = 6e8;// Gb
	// conformational constants
	matParams->ve[13] = 0.1553;// alpha_c
	matParams->ve[14] = 0.001;// eta_c
	matParams->ve[15] = 1.8098e17;// Ns_c
	matParams->ve[16] = 1.38e-17;// boltzmann constant kB
	// slippage
	matParams->ve[17] = 100;// lambdaCrit
	matParams->ve[18] = 383.15;// Ts 
	matParams->ve[19] = 0.359e6;// gamma0_ref = 0.653
	matParams->ve[20] = 7307.8;// Cs 10612
	matParams->ve[21] = 152.95;// Tinf 95.48
	matParams->ve[22] = 0.1565;// C1
	matParams->ve[23] = 39.937;// C2
	matParams->ve[24] = 0.9878;// beta
	matParams->ve[25] = 0.33;// poissons ratio


	// crit lambda properties
	matParams->ve[26] = -0.0111; // C1
	matParams->ve[27] = 3.627; // C2
	matParams->ve[28] = 0.9856; // BETA
	matParams->ve[29] = -0.0356; // k
	matParams->ve[30] = 15.393; // b 



	double t_n_1 = 0;
	double t_n = 0;

	int writeFreq = 10;


	int n = 0;


	state_variables ** stateOld = new_material_state(&TEMPERATURE, 1, 
		BUCKLEY_MATERIAL , 
		PLASTIC_MATERIAL , 
		DIM, IS_AXI);
	state_variables ** stateNew = new_material_state(&TEMPERATURE, 1, 
		BUCKLEY_MATERIAL , 
		PLASTIC_MATERIAL , 
		DIM, IS_AXI);


	VEC * stressVoigt = v_get(5);

	MAT * P = m_get(3,3);
	MAT * inter = m_get(3,3);
	MAT * sigma = m_get(3,3);
	FILE * fp;
	int IS_AXI = 0;

	fp = fopen("Stress_11_nominal.txt","w");
	double maxStrain = 0;

	fclose(fp);
	MAT * invF = m_get(3,3);
	struct timeval start, end;
	gettimeofday(&start, NULL);

	double sigplot = 0 ;
	double strainplot = 0;

	while (( t_n < TMAX) && (maxStrain < PEAK_STRAIN) )
	//while ( n < 1)
	{

		// update time
		t_n_1  = t_n +  DT;

		if ( strcmp (DEFORMATION_MODE, "BIAXIAL") == 0 )
		{
			stateNew[0]->F->me[0][0] = 1.00+t_n_1*SR;
			stateNew[0]->F->me[1][1] = 1.00+t_n_1*SR;
			stateNew[0]->F->me[2][2] = (1.00/pow(1.00+t_n_1*SR,2));

		}else if ( strcmp (DEFORMATION_MODE, "SIMPLE_SHEAR") == 0 ) {

			stateNew[0]->F->me[0][1] = SR*t_n_1; 
			stateNew[0]->F->me[1][0] = 0; 
		}



		/* ------------------------------------------*/
		/* -----------------Deformation--------------*/
		/* ------------------------------------------*/

		/*  Find Rate of deformation tensors at n+h*/
		velocity_grad(stateNew[0],stateOld[0],DT,0.500);

		// Find inverse deformation gradient at n+1
		m_inverse_small(stateNew[0]->F, stateNew[0]->invF);

		// Find Jacobian at n+1
		stateNew[0]->Jacobian = determinant(stateNew[0]->F);

		//------------------------------------------//
		//        Update Polar Decomposition        //
		//------------------------------------------//
		update_Polar_Decomposition(stateNew[0], stateOld[0], DT);

		// remove rotation from D to get d
		un_rotate_tensor(stateNew[0]->D, stateNew[0]->R, 
			stateNew[0]->m_temp1, stateNew[0]->d);


		un_rotate_tensor(stateNew[0]->Omega, stateNew[0]->R, 
			stateNew[0]->m_temp1, stateNew[0]->omega);


		// Find stress
		buckleyStress(stateNew[0], stateOld[0], matParams,DT);




		if ( strcmp (DEFORMATION_MODE, "BIAXIAL") == 0 )
		{
			strainplot = stateNew[0]->F->me[0][0]-1;
			maxStrain = stateNew[0]->F->me[0][0]-1;
			sigplot = (stateNew[0]->Sc->me[0][0] - stateNew[0]->Sc->me[2][2])
			+(stateNew[0]->Sb->me[0][0] - stateNew[0]->Sb->me[2][2]);
			sigplot = sigplot/pow(10,6);

		}else if ( strcmp (DEFORMATION_MODE, "SIMPLE_SHEAR") == 0 ) {
			maxStrain = stateNew[0]->F->me[0][1];
			strainplot = stateNew[0]->F->me[0][1];
			sigplot = stateNew[0]->Sc->me[0][1] + stateNew[0]->Sb->me[0][1] ;//+ stateNew[0]->Sc->me[0][1] ;
			sigplot = sigplot/pow(10,6);
		}
		// write stress to file 
		if ( n % writeFreq == 0)
		{
			fp = fopen("Stress_11_nominal.txt","a");
			//fprintf(fp,"%lf,%lf\n",t_n_1, sig12);
			fprintf(fp,"%lf,%lf\n",strainplot, sigplot);

			fclose(fp);
		}



		//-------------------------------------------------//
		//                  State updates                  //
		//-------------------------------------------------//

		// n+1 --->>>> n 
		m_copy(stateNew[0]->F,stateOld[0]->F);
		m_copy(stateNew[0]->V,stateOld[0]->V);
		m_copy(stateNew[0]->R,stateOld[0]->R);
		m_copy(stateNew[0]->invF,stateOld[0]->invF);			
		stateOld[0]->Jacobian = stateNew[0]->Jacobian;
		stateOld[0]->div_v = stateNew[0]->div_v;



		++n;
		t_n = t_n_1;
	}








	gettimeofday(&end, NULL);
	double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
		end.tv_usec - start.tv_usec) / 1.e6;

	// get time taken to run
	printf("buckley took %lf seconds to run\n", delta);










	// printf("omega =");
	// m_foutput(stdout, stateNew[0]->Omega);
	// printf("V =");
	// m_foutput(stdout, stateNew[0]->V);
	// printf("D =");
	// m_foutput(stdout, stateNew[0]->D);
	// printf("R =");

	// m_foutput(stdout, stateNew[0]->R);


	// printf("det V  = %lf\n", determinant(stateNew[0]->V));
	// printf("det F  = %lf\n", determinant(stateNew[0]->F));

	// m_foutput(stdout, stateNew[0]->F);

	// printf("trace of D = %lf", stateNew[0]->Dbar->me[0][0] + stateNew[0]->Dbar->me[1][1] + stateNew[0]->Dbar->me[2][2]);
	// gettimeofday(&end, NULL);
	// double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
	// 	end.tv_usec - start.tv_usec) / 1.e6;



	//test some eigen value computations

	MAT * test = m_get(3,3);
	test->me[0][0] = 3;
	test->me[0][1] = 2;
	test->me[0][2] = 4;

	test->me[1][0] = 2;
	test->me[1][1] = 0;
	test->me[1][2] = 2;

	test->me[2][0] = 4;
	test->me[2][1] = 2;
	test->me[2][2] = 3;


	MAT * Q = m_get(3,3); 

	VEC * eigVals = v_get(3);

	
	dsyevq3(test->me, Q->me, eigVals->ve);

	MAT * temp = m_get(3,3);
	MAT * temp_1 = m_get(3,3);


	temp->me[0][0] = eigVals->ve[0];
	temp->me[1][1] = eigVals->ve[1];
	temp->me[2][2] = eigVals->ve[2];

	m_mlt(Q,temp,temp_1);
	mmtr_mlt(temp_1,Q,temp);


	m_foutput(stdout, temp);
	v_foutput(stdout, eigVals);

	m_foutput(stdout, Q);

	// MAT * m_temp = m_get(3,3);
	// MAT * m_temp1 = m_get(3,3);


	// MAT * D = m_get(3,3);

	// D->me[0][0] = eigVals->ve[0];
	// D->me[1][1] = eigVals->ve[1];
	// D->me[2][2] = eigVals->ve[2];

	// m_mlt(test,Q,m_temp1);
	// m_mlt(Q,D,m_temp);

	// m_sub(m_temp1,m_temp,m_temp);

	// m_foutput(stdout, m_temp);




	PERM * order = px_get(3);
	v_sort(eigVals, order);
	px_free(order);


	exit(0);



}