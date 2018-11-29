#include "weight_function.h"



#include <math.h>

int weight_function (VEC * weights, double  * xS, double dI, char * type,  int compute, int dim){


	double r = 0;
	for ( int i = 0 ; i < dim ; i++){
		r += pow(xS[i],2);
	}

	r = sqrt(r)/dI;

	if ( r > 1)
	{
		printf("ERROR R > 1 \n");
		printf("r = %lf\n",r);
		printf("xS[0] = %lf xS[1] = %lf \n ", xS[0], xS[1]);
	}

	// check input matrix is of the right dimensions 
	if (weights == VNULL){
		if (compute == 1){
			weights = v_get(1);
		}else if ( compute == 2){
			weights = v_get(1+dim);
		}else if ( compute == 3){
			weights = v_get(1+dim + dim*dim);
		}else{
			fprintf(stderr,"compute variable not valid");
		}
	}else
	{
		if (compute == 1){

			if (weights->max_dim != 1)
			{
				v_resize(weights, 1);
			}
		}else if ( compute == 2){
			if (weights->max_dim != 1+dim)
			{
				v_resize(weights, 1+dim);
			}
		}else if ( compute == 3){
			if (weights->max_dim != 1+dim+dim*dim)
			{
				v_resize(weights, 1+dim+dim*dim);
			}
		}else{
			fprintf(stderr,"compute variable not valid");
		}
	}


	if (strcmp(type, "cubic")  == 0){
		VEC * w_arr = v_get(3);

		// find w, w_r, w_rr 
		cubic_spline(w_arr,r);

		weights->ve[0] = w_arr->ve[0];

		if (( compute == 2) || ( compute == 3)){
			for ( int i = 0 ; i < dim ; i++)
			{
				if ( r != 0 ){
				double drdi = (xS[i]/(r*dI*dI));
				weights->ve[1+i] = drdi*w_arr->ve[1];
				}
				else{
				weights->ve[1+i] = 0;
				}
			}

		}
		if ( compute == 3){
			printf("use higher order basis to find 2nd derivative \n");
		}


	}else if ( strcmp(type,"quartic") == 0)
	{
		VEC * w_arr = v_get(3);

		// find w, w_r, w_rr 
		quartic_spline(w_arr,r);

		weights->ve[0] = w_arr->ve[0];

		if (( compute == 2) || ( compute == 3)){
			for ( int i = 0 ; i < dim ; i++)
			{
				if ( r != 0 ){
				double drdi = (xS[i]/(r*dI*dI));
				weights->ve[1+i] = drdi*w_arr->ve[1];
				}
				else{
				weights->ve[1+i] = 0;
				}
			}

		}

		if ( compute == 3){

			double w_r = w_arr->ve[1];
			double w_rr = w_arr->ve[2];
			int kron = 0;

			for ( int i = 0 ; i < dim ; i++)
			{
				for ( int j = 0 ; j < dim ; j++)
				{
					int indx = j + i*dim;
					if ( i == j ){
						kron = 1;
					}else{
						kron = 0;
					}
					if ( r != 0 )
					{
						weights->ve[(1+dim)+indx] = (w_rr - w_r/r) * ( xS[i]*xS[j] )/(r*r*pow(dI,4))  + w_r * (double)(kron/(r*dI*dI)) ;	
	
					}else{
						weights->ve[(1+dim)+indx] = -12*kron/(dI*dI);	
	
					}
				}
			}
		}

		v_free(w_arr);

		
	}else {
		fprintf(stderr,"weight function not set \n ");
	}



	return 0;
}