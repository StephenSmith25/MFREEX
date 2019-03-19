#include "ShapeFunction/weight_function_materialpoint.h"
#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )


#include <math.h>

static inline int cubic_spline_a(double w_arr[3], double r)
{


    if (r <= 0.50 ){
        
        w_arr[0] = 2.0/3 - 4*r*r + 4*r*r*r;
        w_arr[1] = -8*r + 12*r*r;
    }
    else if ((r > 0.50) & (r <= 1.00)){
        w_arr[0] = 4.0/3 - 4*r + 4*r*r - (4.0/3)*r*r*r;
        w_arr[1] = -4 + 8*r - 4*r*r;
    }
    else{        
        w_arr[0] = 0 ;
        w_arr[1] = 0 ;
        
    }
        
    return 0 ;
}

static inline int quartic_spline_a(double w_arr[3], double r)
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

	w_arr[0] = w;
	w_arr[1] = w_r;
	w_arr[2] = w_rr;


	return 0; 
}

int weight_function_materialpoint(VEC * weights, MATERIAL_POINT * MP, double  * xS, int dim)
{





	char * weight_type = "cubic";

	int m,n;

	double d_I_s = 0;
	double r_norm_s = 0;
	double r_norm;
	double distance = 0;
	double inter[dim] ;

	double dI = MP->r_cutoff;


	/* Find r and dr/dk */
	double drdk[dim];
	switch (dim)
	{

		case(2):
		{
			r_norm = sqrt(xS[0]*xS[0] + xS[1]*xS[1])/dI;

			drdk[0] = (1.00/(r_norm*dI *dI)) * (xS[0]) ;
			drdk[1] = (1.00/(r_norm*dI*dI)) * (xS[1]) ;



			break;
		}
		case(3):
		{
			// r_norm_s = xS[0]*(MP->invMI->me[0][0]*xS[0] + MP->invMI->me[0][1]*xS[1] + MP->invMI->me[0][2]*xS[2] ) + 
			// xS[1]*(MP->invMI->me[1][0]*xS[0] + MP->invMI->me[1][1]*xS[1] + MP->invMI->me[1][2]*xS[2]) +
			// xS[2]*(MP->invMI->me[2][0]*xS[0] + MP->invMI->me[2][1]*xS[1] + MP->invMI->me[2][2]*xS[2]) ;
			// r_norm = sqrt(r_norm_s);
			// // drdk
			// drdk[0] = (1.00/(2*r_norm)) * (2*MP->invMI->me[0][0]*xS[0] + MP->invMI->me[0][1]*xS[1] + MP->invMI->me[0][2]*xS[2]) ;
			// drdk[1] = (1.00/(2*r_norm)) * (MP->invMI->me[1][0]*xS[0] + 2*MP->invMI->me[1][1]*xS[1] + MP->invMI->me[1][2]*xS[2]); 
			// drdk[2] = (1.00/(2*r_norm)) * (MP->invMI->me[2][0]*xS[0] + MP->invMI->me[2][1]*xS[1] + 2*MP->invMI->me[2][2]*xS[2]);

		}
	}





	/* Evaluate weight function */
	if (strcmp(weight_type, "cubic")  == 0){
		double w_arr[3] ;

		// find w, w_r, w_rr 
		cubic_spline_a(w_arr,r_norm);

		weights->ve[0] = w_arr[0];

		for ( int i = 0 ; i < dim; i++)
		{
			if ( r_norm != 0 ){
				double drdi = drdk[i];
				weights->ve[1+i] = drdi*   w_arr[1];
			}
			else{
				weights->ve[1+i] = 0;
			}
		}



	}else if ( strcmp(weight_type,"quartic") == 0)
	{
		double w_arr[3] ;

		// find w, w_r, w_rr 
		quartic_spline_a(w_arr,r_norm);

		weights->ve[0] = w_arr[0];

		for ( int i = 0 ; i < dim ; i++)
		{
			if ( r_norm != 0 ){
				double drdi = drdk[i];
				weights->ve[1+i] = drdi*w_arr[1];
			}
			else{
				weights->ve[1+i] = 0;
			}
		}

	}else{
		printf("weight function not set \n");
	}

	return 0;
}