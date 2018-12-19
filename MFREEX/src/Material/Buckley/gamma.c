/*
 * =====================================================================================
 *
 *       Filename:  gamma.c
 *
 *    Description:  Finds the viscosity of the dashpot used in the conformational part of
 *    		    the Buckley model. 
 *
 *        Version:  1.0
 *        Created:  10/05/18 17:56:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  STEPHEN SMITH 	 
 *   Organization:  QUEEN'S UNIVERSITY BELFAST 
 *
 * =====================================================================================
 */

#include "Material/Buckley/gamma.h"

double gammaV(state_Buckley * state, double maxLambdaN,double critLambda,VEC * para){


	double starT = para->ve[18];
	double refGamma = para->ve[19];
	double Cs = para->ve[20];
	double Tinf = para->ve[21];
	double shiftTSlope = para->ve[22];
	double shiftTIntercept = para->ve[23];
	double expFactor = para->ve[24];
	double gamma0 = 0;
	VEC * lambda = state->eigValVBar;
	VEC * eigD = state->eigValDBar;



	double temperature = state->temperature;


	// initialise gamma 
	double gamma_n_1 = 0; 

	// find theta 
	VEC * strainRate = v_get(3);
	int index = 0;
	double theta = 0;
	// find the strain rate
	for ( int i = 0 ; i < 3 ; i++){
		strainRate->ve[i] = eigD->ve[i]*lambda->ve[i] ;
	}


	double maxSr = v_max(strainRate,&index);



	if ( maxSr == 0){
		maxSr = 0.0000001; 
	}

	// find theta, the ratio of in plane natural strain rate
	if ( eigD->ve[1] < 1e-6 ){
		theta = 0;
	}else{
		theta = eigD->ve[2]/eigD->ve[1];
	}

	// find xi 
	double xi = ( 2 * theta + 1)/(theta +2);
	if ( xi > 1){
		xi = 1;
	}else if( xi < 0 ) {
		xi = 0.0001;
	}else{
		// do nothing
	}

	double log2sr = log(maxSr)/log(2);
	if ( maxSr < 0.001){
	log2sr = 0;
	}


	double shiftTemperature = temperature * pow(10, ( ( shiftTSlope*log2sr ) / ( shiftTIntercept + log2sr) )* pow(expFactor,2-2*xi));


	gamma0 = exp( Cs/(shiftTemperature - Tinf) - Cs/(starT - Tinf));
	gamma0 = gamma0*refGamma;
	if ( maxLambdaN >= critLambda ){
		gamma_n_1 = 1e50;

	}else{
		gamma_n_1 = gamma0/ ( 1.000 - (maxLambdaN/critLambda)) ; 
	}
	
	v_free(strainRate);
	
	
	gamma_n_1 = gamma_n_1 ;
	
	return gamma_n_1;
	


}

