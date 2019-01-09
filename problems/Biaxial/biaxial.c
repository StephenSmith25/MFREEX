

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "Material/Buckley/new_Buckley_State.h"
#include "Material/Buckley/buckleyStress.h"
#include "Material/material.h"
#include <sys/time.h>

int main(void)
{

	double temperature = 85;
	double sr = 4;
	double peakStrain = 2.5;

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



	double dt = 1e-4;

	double tmax = 5;
	double t_n_1 = 0;
	double t_n = 0;

	int writeFreq = 10;


	int n = 0;
	state_Buckley ** stateOld = new_Buckley_State(1, &temperature, is_axi, dim);
	state_Buckley ** stateNew = new_Buckley_State(1, &temperature, is_axi, dim);



	FILE * fp;

	fp = fopen("Stress_11_nominal.txt","w");
	double maxStrain = 0;

	fclose(fp);

	struct timeval start, end;
	gettimeofday(&start, NULL);


	while (( t_n < tmax) && (maxStrain < peakStrain))
	{
		t_n_1  = t_n +  dt;


		stateNew[0]->F->me[0][0] = 1.00+t_n_1*sr;
		stateNew[0]->F->me[1][1] = 1.00+t_n_1*sr;
		stateNew[0]->F->me[2][2] = 1.00/pow(1.00+t_n_1*sr,2);

		buckleyStress(stateNew[0], stateOld[0], matParams,critLambdaParams,dt);



		double sig11 = stateNew[0]->Sc->me[0][0] - stateNew[0]->Sc->me[2][2]+stateNew[0]->Sb->me[0][0] - stateNew[0]->Sb->me[2][2];
		sig11 = sig11/pow(10,6);
		maxStrain = stateNew[0]->F->me[0][0]-1;
		// write stress to file 
		if ( n % writeFreq == 0)
		{
			fp = fopen("Stress_11_nominal.txt","a");
			fprintf(fp,"%lf,%lf\n",stateNew[0]->F->me[0][0]-1, sig11);
			fclose(fp);
		}
		++n;
		t_n = t_n_1;
	}


	gettimeofday(&end, NULL);
	double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6;

	// get time taken to run
	printf("buckley took %lf seconds to run\n", delta);










	exit(0);



}