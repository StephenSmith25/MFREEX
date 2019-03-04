#include "weight_function.h"
#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )


#include <math.h>


int weight_function (VEC * weights, double  * xS, int I, meshfreeDomain * mfree,  int compute){

	int dim = mfree->dim;

	char * shape = mfree->kernel_shape;
	char * weight_type = mfree->weight_function;
	enum SUPPORT_TYPE support = mfree->kernel_support;



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

	double dI = mfree->di->ve[I];

	double r = 0;
	for ( int i = 0 ; i < dim ; i++){
		r += pow(xS[i],2);
	}

	r = sqrt(r)/dI;

	if (strcmp(weight_type, "cubic")  == 0){
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



	}else if ( strcmp(weight_type,"quartic") == 0)
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

	double dI_i[2] = {mfree->di_tensor->me[I][0],mfree->di_tensor->me[I][1]};
	double r_j[3];

	for ( int i = 0 ; i < dim ; i++)
	{
		r_j[i] = fabs(xS[i])/dI_i[i];
	}


	// else if tensor product kernels
	
	if ( strcmp(weight_type,"cubic") == 0)
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

	double d_I_s = 0;
	MAT * MI = mfree->MI[I];


	double M11 = MI->me[0][0];
	double M12 = MI->me[0][1];
	double M21 = MI->me[1][0];
	double M22 = MI->me[1][1];



	double r_norm_s = xS[0]*(M11*xS[0] + M12 * xS[1]) + xS[1]*(M21*xS[0] + M22 * xS[1]);
	double r_norm = sqrt(r_norm_s);

	double drdk[2];

	drdk[0] = (1.00/(2*r_norm) * (2*M11*xS[0] + M12*xS[1]))  ; 
	drdk[1] = (1.00/(2*r_norm) * (2*M22*xS[1] + M21*xS[0])) ; 


	if (strcmp(weight_type, "cubic")  == 0){
		VEC * w_arr = v_get(3);

		// find w, w_r, w_rr 
		cubic_spline(w_arr,r_norm);

		weights->ve[0] = w_arr->ve[0];

		if (( compute == 2) || ( compute == 3)){
			for ( int i = 0 ; i < 2; i++)
			{
				if ( r_norm != 0 ){
				double drdi = drdk[i];
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



	}else if ( strcmp(weight_type,"quartic") == 0)
	{
		VEC * w_arr = v_get(3);

		// find w, w_r, w_rr 
		quartic_spline(w_arr,r_norm);

		weights->ve[0] = w_arr->ve[0];

		if (( compute == 2) || ( compute == 3)){
			for ( int i = 0 ; i < dim ; i++)
			{
				if ( r_norm != 0 ){
				double drdi = drdk[i];
				weights->ve[1+i] = drdi*w_arr->ve[1];
				}
				else{
				weights->ve[1+i] = 0;
				}
			}

		}

		// if ( compute == 3){

		// 	double w_r = w_arr->ve[1];
		// 	double w_rr = w_arr->ve[2];
		// 	int kron = 0;

		// 	for ( int i = 0 ; i < dim ; i++)
		// 	{
		// 		for ( int j = 0 ; j < dim ; j++)
		// 		{
		// 			int indx = j + i*dim;
		// 			if ( i == j ){
		// 				kron = 1;
		// 			}else{
		// 				kron = 0;
		// 			}
		// 			if ( r_norm != 0 )
		// 			{
		// 				weights->ve[(1+dim)+indx] = (w_rr - w_r/r_norm) * ( xS[i]*xS[j] )/(r_norm*r_norm*pow(dI,4))  + w_r * (double)(kron/(r_norm*dI*dI)) ;	
	
		// 			}else{
		// 				weights->ve[(1+dim)+indx] = -12*kron/(dI*dI);	
	
		// 			}
		// 		}
		// 	}
		// }

		v_free(w_arr);



	}else {
		fprintf(stderr,"weight function not set \n ");
	}
	

		break;
	}
	}


	return 0;
}