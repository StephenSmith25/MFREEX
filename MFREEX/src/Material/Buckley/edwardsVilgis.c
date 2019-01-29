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

int edwardsVilgis(VEC * stress, VEC * lambda, VEC * para, double Jacobian, double temperature){


	double alpha = para->ve[13];
	double Ns = para->ve[15];
	double eta = para->ve[14];
	double kB = para->ve[16];


	double sigCL = 0;
	double sigSL = 0; 



	double I1 = lambda->ve[0] + lambda->ve[1] + lambda->ve[2];


	double ta = Ns*kB*temperature;
	double tb = 0;
	double tc = 1 + eta;
	double td = pow(alpha,2);
	double te = 1 - td;
	double tf = 1- td * I1;
	double tg = 1 - 2*td ;
	double th = 0;


	for ( int k = 0 ; k < 3 ; k++){
		th = th + (lambda->ve[k])/(1+eta*lambda->ve[k]);
	}




	for ( int k = 0 ; k < 3 ; k++){



		// sigCL part 
		double tempa = tg/tf;
		double tempb = td * te * I1/pow(tf,2);
		double tempc = lambda->ve[k] * tb;
		sigCL = tempc * ( tempa + tempb);


		// sigSL 
		tempa = lambda->ve[k] * ta;
		tempb = tc*te;
		tempc = 1 + eta*lambda->ve[k];
		double tempd = tempb/(tf*pow(tempc,2));
		double tempe = tempb * td * th / pow(tf,2);
		double tempf = eta/tempc;
		double tempg = td/tf;
		sigSL = tempa * ( tempd + tempe + tempf - tempg);
		stress->ve[k] = (1/Jacobian)*(sigCL + sigSL) * 1e3;


	}

	double I1s = (1.000/3.000)*(stress->ve[0] + stress->ve[1] + stress->ve[2]);
	for ( int k = 0 ; k < 3 ; k++){
		stress->ve[k] = stress->ve[k] - I1s;
	}



	return 0 ;
}
