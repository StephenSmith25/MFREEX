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


int edwardsVilgis(VEC * stress, VEC * lambda_in, VEC * para, double Jacobian, double temperature){



	// Get Edwards-vilgis consants from parameter inputs
	double alpha = para->ve[13];
	double Ns = para->ve[15];
	double eta = para->ve[14];
	double kB = para->ve[16];


	// Convert log stretch into stretch
	double lambda1 = exp(lambda_in->ve[0]);
	double lambda2 = exp(lambda_in->ve[1]);
	double lambda3 = exp(lambda_in->ve[2]);

	// find square of stretches
	lambda1 = lambda1*lambda1;
	lambda2 = lambda2*lambda2;
	lambda3 = lambda3*lambda3;

	// create array of stretches
	double lambda[3] = {lambda1,lambda2,lambda3};



	// Constants 
	double D = Ns * kB * temperature*0.5;
	double E = (1+eta)*(1-alpha*alpha);
	double zeta = 1 - alpha*alpha * (lambda1  + lambda2  + lambda3);
	double xi = lambda1/(1+eta*lambda1) + 	lambda2/(1+eta*lambda2) +lambda3/(1+eta*lambda3);
	double phi = log(1 + eta*lambda1) + log(1+eta*lambda2) + log(1 + eta*lambda3);


	// First invariant of stretches
	double I1 = lambda1 + lambda2 + lambda3;

	// Form each principle component of stress
	for ( int k = 0 ; k < 3 ; k++){


		stress->ve[k] = (1000/Jacobian)*(  
		2*D*alpha*alpha*lambda[k] * (E*xi/(zeta*zeta) - 1.00/zeta)
		+ D*E * 2*lambda[k] / (zeta*pow((1+eta*lambda[k]),2))
		+ D * eta*lambda[k]/(1+eta*lambda[k])   
		);


	}

	// Subtract hydrostatic pressure term
	double I1s = (1.000/3.000)*(stress->ve[0] + stress->ve[1] + stress->ve[2]);
	for ( int k = 0 ; k < 3 ; k++){
		stress->ve[k] = stress->ve[k] - I1s;
	}

	return 0 ;
}

// 	double alpha = para->ve[13];
// 	double Ns = para->ve[15];
// 	double eta = para->ve[14];
// 	double kB = para->ve[16];

// 	double lambda1 = exp(lambda_in->ve[0]);
// 	double lambda2 = exp(lambda_in->ve[1]);
// 	double lambda3 = exp(lambda_in->ve[2]);


// 	lambda1 = lambda1*lambda1;
// 	lambda2 = lambda2*lambda2;
// 	lambda3 = lambda3*lambda3;



// 	double lambda[3] = {lambda1,lambda2,lambda3};

// 	double sigCL = 0;
// 	double sigSL = 0; 



// 	double I1 = lambda[0] + lambda[1] + lambda[2];


// 	double ta = Ns*kB*temperature;
// 	double tb = 0;
// 	double tc = 1 + eta;
// 	double td = pow(alpha,2);
// 	double te = 1 - td;
// 	double tf = 1- td * I1;
// 	double tg = 1 - 2*td ;
// 	double th = 0;




// 	for ( int k = 0 ; k < 3 ; k++){
// 		th = th + (lambda[k])/(1+eta*lambda[k]);
// 	}



// 	for ( int k = 0 ; k < 3 ; k++){



// 		// sigCL part 
// 		double tempa = tg/tf;
// 		double tempb = td * te * I1/pow(tf,2);
// 		double tempc = lambda[k] * tb;
// 		sigCL = tempc * ( tempa + tempb);



// 		// sigSL 
// 		tempa = lambda[k] * ta;

// 		tempb = tc*te;
// 		tempc = 1 + eta*lambda[k];
// 		double tempd = tempb/(tf*pow(tempc,2));
// 		double tempe = tempb * td * th / pow(tf,2);
// 		double tempf = eta/tempc;


// 		double tempg = td/tf;
// 		if ( lambda[k] == 1.00 )
// 		{
// 			tempa = 0;
// 		}

// 		sigSL = tempa * ( tempd + tempe + tempf - tempg);

// 		stress->ve[k] = (1/Jacobian)*(sigCL + sigSL) * 1e3;


// 	}

// 	double I1s = (1.000/3.000)*(stress->ve[0] + stress->ve[1] + stress->ve[2]);

// 	for ( int k = 0 ; k < 3 ; k++){
// 		stress->ve[k] = stress->ve[k] - I1s;
// 	}




// 	return 0 ;
// }