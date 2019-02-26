#include "weight_function.h"
#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )


#include <math.h>


int weight_function (VEC * weights, double  * xS, double * dI_i, char * type, enum SUPPORT_TYPE support,  int compute, int dim){

	double dI = dI_i[0];

	double r = 0;
	for ( int i = 0 ; i < dim ; i++){
		r += pow(xS[i],2);
	}

	r = sqrt(r)/dI;

	// if tensor product is used
	double r_j[3];

	for ( int i = 0 ; i < dim ; i++)
	{
		r_j[i] = fabs(xS[i])/dI_i[i];
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






	switch (support)
	{

	case (RADIAL):
	{


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
		v_free(w_arr);



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
	

		break;
	}


	// using rectangular kernels ( tensor products)

	case (RECTANGULAR):
	{


	// else if tensor product kernels
	
	if ( strcmp(type,"cubic") == 0)
	{
		VEC * w_arr = v_get(3);
		weights->ve[0] = 1;

		double w_j[3];

		for ( int i = 0 ; i < dim ; i++)
		{
			cubic_spline(w_arr,r_j[i]);
			weights->ve[0] = weights->ve[0]*w_arr->ve[0];
			w_j[i] = w_arr->ve[0];

		}


		if ( compute == 2){
			for ( int i = 0 ; i < dim ; i++)
			{

				printf("getting first derivative\n \n");
				double r = r_j[i];
				double wjwk  = 1;
				for ( int j = 0 ; j < dim ; j++)
				{
					if ( j != i)
					{
						wjwk = wjwk*w_j[j];
					}
				}


				if ( r != 0 ){
					double drdi = sign(xS[i])/(dI_i[i]);

					printf("drdi = %lf \n", drdi);
					weights->ve[1+i] = drdi*w_arr->ve[1]*wjwk;
				}
				else{
					weights->ve[1+i] = 0;
				}
			}

		}
		// Find w_x, and w_y


		v_free(w_arr);



	}else {
		fprintf(stderr,"weight function not set \n ");
	}



		break;
	}
	case (ELLIPTICAL):
	{
		fprintf(stderr,"need to implement this");
	}
	}


	return 0;
}