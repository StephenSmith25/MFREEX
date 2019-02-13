/*
 * =====================================================================================
 *
 *       Filename:  edwardsVilgis.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/05/18 09:30:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#include "Material/Buckley/edwardsVilgis.h"

static int call_count = 0;

int edwardsVilgis(VEC * stress, VEC * lambda_in, VEC * para, double Jacobian, double temperature){


	double alpha = para->ve[13];
	double Ns = para->ve[15];
	double eta = para->ve[14];
	double kB = para->ve[16];

	double lambda1 = exp(lambda_in->ve[0]);
	double lambda2 = exp(lambda_in->ve[1]);
	double lambda3 = exp(lambda_in->ve[2]);


	lambda1 = lambda1*lambda1;
	lambda2 = lambda2*lambda2;
	lambda3 = lambda3*lambda3;



	double lambda[3] = {lambda1,lambda2,lambda3};


	double D = Ns * kB * temperature*0.5;
	double E = (1+eta)*(1-alpha*alpha);


	double zeta = 1 - alpha*alpha * (lambda1  + lambda2  + lambda3);
	double xi = lambda1/(1+eta*lambda1) + 	lambda2/(1+eta*lambda2) +lambda3/(1+eta*lambda3);
	double phi = log(1 + eta*lambda1) + log(1+eta*lambda2) + log(1 + eta*lambda3);



	double I1 = lambda1 + lambda2 + lambda3;


	for ( int k = 0 ; k < 3 ; k++){


		stress->ve[k] = (1000/Jacobian)*(  
		2*D*alpha*alpha*lambda[k] * (E*xi/(zeta*zeta) - 1.00/zeta)
		+ D*E * 2*lambda[k] / (zeta*pow((1+eta*lambda[k]),2))
		+ D * eta*lambda[k]/(1+eta*lambda[k])   
		);


	}


	double I1s = (1.000/3.000)*(stress->ve[0] + stress->ve[1] + stress->ve[2]);
	for ( int k = 0 ; k < 3 ; k++){
		stress->ve[k] = stress->ve[k] - I1s;
	}
	if ( call_count == 0)
	{
		v_foutput(stdout, stress);
	}	
	++call_count;


	return 0 ;
}
