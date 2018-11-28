#include "quartic_spline.h"


int quartic_spline(VEC * w_arr, double r)
{
	double w = 0;
	double w_r = 0;
	double w_rr = 0;

	if ( r <= 1 )
	{
		w = 1 - 6 * pow(r,2) + 8*pow(r,3) - 3*pow(r,4);
		w_r = -12*r + 24*pow(r,2) - 12 * pow(r,3);
		w_rr = -12 + 48*r - 36*pow(r,2);
 
	}else{


	}

	w_arr->ve[0] = w;
	w_arr->ve[1] = w_r;
	w_arr->ve[2] = w_rr;


	return 0; 
}