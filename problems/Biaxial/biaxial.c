

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "Material/Buckley/new_Buckley_State.h"
#include "Material/Buckley/buckleyStress.h"
#include "Material/material.h"
#include <sys/time.h>
#include "Material/Hyperelastic/yeoh.h"
int main(void)
{

	double temperature = 105;
	double sr = 16;
	double peakStrain = 3;

	double dim = 3;
	double is_axi = 0;

	// material parameters 

	/*  Material struct */
	VEC * matParams = v_get(26);
	matParams->ve[0] = 2.814e-3; // VS
	matParams->ve[1] = 0.526e-3; // VP
	matParams->ve[2] = 1.81e6; // mu*_0
	matParams->ve[3] = 328.76; // Tinf
	matParams->ve[4] = 358.15; // T*
	matParams->ve[5] = matParams->ve[4]; // Tf*
	matParams->ve[6] = (67.47); // Cv
	matParams->ve[7] = 1.23e5; // H0
	matParams->ve[8] = 8.314; // R
	matParams->ve[9] = 1.8e9; // Kb
	matParams->ve[10] = 1e7;// Gb
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
	VEC * critLambdaParams = v_get(6);
	critLambdaParams->ve[0] = -0.0111; // C1
	critLambdaParams->ve[1] = 3.627; // C2
	critLambdaParams->ve[2] = 0.9856; // BETA
	critLambdaParams->ve[3] = -0.0356; // k
	critLambdaParams->ve[4] = 15.393; // b 


	VEC * materialParameters = v_get(4);

	materialParameters->ve[0] = 3.0208;
	materialParameters->ve[1] =-0.1478;
	materialParameters->ve[2] = 0.0042;
	materialParameters->ve[3] =100;
	double dt = 1e-5;

	double tmax = 10;
	double t_n_1 = 0;
	double t_n = 0;

	int writeFreq = 10;


	int n = 0;
	state_Buckley ** stateOld = new_Buckley_State(1, &temperature, is_axi, dim);
	state_Buckley ** stateNew = new_Buckley_State(1, &temperature, is_axi, dim);
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

	while (( t_n < tmax) && (maxStrain < peakStrain) )
	{
		t_n_1  = t_n +  dt;


	
		stateNew[0]->F->me[0][0] = 1.00+t_n_1*sr;
		stateNew[0]->F->me[1][1] = 1.00+t_n_1*sr;
		stateNew[0]->F->me[2][2] = (1.00/pow(1.00+t_n_1*sr,2));


		// stateNew[0]->F->me[0][1] = sr*t_n_1; 
		// stateNew[0]->F->me[1][0] = 0; 
	
		

		buckleyStress(stateNew[0], stateOld[0], matParams,critLambdaParams,dt,0, IS_AXI);



		// m_inverse(stateNew[0]->F, invF);

		// yeoh(stressVoigt,stateNew[0]->F,materialParameters);

		// P->me[0][0] = stressVoigt->ve[0];
		// P->me[1][1] = stressVoigt->ve[1];
		// P->me[0][1] = stressVoigt->ve[2];
		// P->me[1][0] = stressVoigt->ve[3];
		// P->me[2][2] = stressVoigt->ve[4];

		// mmtr_mlt(P, stateNew[0]->F, sigma);







		double sig11 = (stateNew[0]->Sc->me[0][0] - stateNew[0]->Sc->me[2][2])
		+(stateNew[0]->Sb->me[0][0] - stateNew[0]->Sb->me[2][2]);
		sig11 = sig11/pow(10,6);
		//sig11 = stateNew[0]->Sc->me[0][0] + stateNew[0]->Sb->me[0][0];
		//sig11 = sig11/pow(10,6);
		double sig12 = 0*stateNew[0]->Sc->me[0][1] + stateNew[0]->Sb->me[0][1] ;//+ stateNew[0]->Sc->me[0][1] ;
		sig12 = sig12/pow(10,6);

		//sig11 = sigma->me[1][1];
		maxStrain = stateNew[0]->F->me[0][0]-1;
		//maxStrain = stateNew[0]->F->me[0][1];
		// write stress to file 
		if ( n % writeFreq == 0)
		{
			fp = fopen("Stress_11_nominal.txt","a");
			//fprintf(fp,"%lf,%lf\n",t_n_1, sig12);
			fprintf(fp,"%lf,%lf\n",stateNew[0]->F->me[0][0]-1, sig11);

			fclose(fp);
		}
		++n;
		t_n = t_n_1;
	}

	printf("omega =");
	m_foutput(stdout, stateNew[0]->Omega);
	printf("V =");
	m_foutput(stdout, stateNew[0]->V);
	printf("D =");
	m_foutput(stdout, stateNew[0]->D);
	printf("R =");

	m_foutput(stdout, stateNew[0]->R);


	printf("det V  = %lf\n", determinant(stateNew[0]->V));
	printf("det F  = %lf\n", determinant(stateNew[0]->F));

	m_foutput(stdout, stateNew[0]->F);

	printf("trace of D = %lf", stateNew[0]->Dbar->me[0][0] + stateNew[0]->Dbar->me[1][1] + stateNew[0]->Dbar->me[2][2]);
	gettimeofday(&end, NULL);
	double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6;

	// get time taken to run
	printf("buckley took %lf seconds to run\n", delta);


	// test some eigen value computations


	MAT * test = m_get(3,3);
	test->me[0][0] = 1;
	test->me[0][1] = 88;
	test->me[0][2] = 150;

	test->me[1][0] = 31.1;
	test->me[1][1] = 1;
	test->me[1][2] = 95.1;

	test->me[2][0] = 11.4;
	test->me[2][1] = 13.1;
	test->me[2][2] = 1;
	MAT * Q = m_get(3,3); 
	VEC * eigVals = v_get(3);

	symmeig(test, Q, eigVals);
	PERM * order = px_get(3);
	v_sort(eigVals, order);
	px_free(order);
	v_foutput(stdout, eigVals);


	exit(0);



}