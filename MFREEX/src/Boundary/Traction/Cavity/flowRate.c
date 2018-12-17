/*
 * =====================================================================================
 *
 *       Filename:  flowRate.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  22/03/18 11:50:24
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include "Boundary/Traction/Cavity/flowRate.h"


double flowRate(double pRatio,double T0, double P0,double R, double Ar, double gamma )
{ 	

	if ( pRatio > 1){
		pRatio =1;
	}
	double massFlow = (Ar*P0/sqrt(R*T0))*pow(pRatio,1.000/gamma)*sqrt(((2*gamma)/(gamma-1))*(1-pow(pRatio,(gamma-1)/gamma))); 

	return massFlow;
}
