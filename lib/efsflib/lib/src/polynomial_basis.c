

#include "polynomial_basis.h"


/* Returns the polynomial basis at the point (double * point) with dim
	order can be linear or quadratic basis
	compute can be 1 (just p(x)) compute = 2 ( dpdx_i) compute = 3 (d2p dx_j dx_i)
	basis terms returned in the columns of the matrix basis
*/

int polynomial_basis(MAT * basis, double * point, int dim, char * order, int compute)
{

	if ( basis != MNULL)
	{
		m_zero(basis);	
	}
	// check basis is right dimensions
	// get size of matricies based on basis used

	int dim_p;

	
	if ( strcmp(order,"linear") == 0 )
	{
		if ( dim == 1){
			dim_p = 2;
		}else if ( dim == 2){
			dim_p = 3;
		}else{
			dim_p = 4;
		}
	}else if ( strcmp(order,"quadratic") == 0)
	{
		if ( dim == 1){
			dim_p = 3;
		}else if ( dim == 2){
			dim_p = 6;
		}else{
			dim_p = 10;
		}

	}else {
		fprintf(stderr,"basis type not specified, or unknown");
	}	

	int dim_c = 1 + dim + dim*dim;

	
	if (strcmp(order, "linear") == 0){


		if ( compute > 2)
		{
			fprintf(stderr, "compute variable too high, use higher order basis\n");
		}
		// check matrix is corrrect size
		if (( basis->m != dim_p ) || (basis->n != dim_c))
		{
			m_resize(basis,dim_p,dim_c);
		}

		if ( dim == 1)
		{
			// p(x)
			basis->me[0][0] = 1;
			basis->me[1][0] = point[0];

			if ( compute == 2){
			// dpdx
				basis->me[1][1] = 1;
			}else {

			}

		}else if ( dim == 2)
		{
			// p(x)
			basis->me[0][0] = 1;
			basis->me[1][0] = point[0];
			basis->me[2][0] = point[1];

			if ( compute == 2){
				// dpdx
				basis->me[1][1] = 1;
				// dpdy
				basis->me[2][2] = 1;
			}else {

			}


		}else if (dim == 3){
			// p(x)
			basis->me[0][0] = 1;
			basis->me[1][0] = point[0];
			basis->me[2][0] = point[1];
			basis->me[3][0] = point[2];

			if ( compute == 2){
				// dpdx
				basis->me[1][1] = 1;
				// dpdy
				basis->me[2][2] = 1;
				// dpdz
				basis->me[3][3] = 1;
			}else {

			}


		}else{
			fprintf(stderr,"dimension not set\n");
		}

	}else if (strcmp(order, "quadratic") == 0){

		// check matrix is corrrect size
		if (( basis->m != dim_p ) || (basis->n != dim_c ))
		{
			m_resize(basis,dim_p,dim_c);
		}

		if ( dim == 1)
		{
			// p(x)
			basis->me[0][0] = 1;
			basis->me[1][0] = point[0];
			basis->me[1][0] = point[0]*point[0];

			if ( compute == 2){
			// dpdx
				basis->me[1][1] = 1;
				basis->me[2][1] = 2*point[0];
			}else if ( compute == 3){
				basis->me[2][2] = 2;
			}else{

			}

		}else if ( dim == 2){
			
			// p(x)
			basis->me[0][0] = 1;
			basis->me[1][0] = point[0];
			basis->me[2][0] = point[1];
			basis->me[3][0] = point[0]*point[0];
			basis->me[4][0] = point[1]*point[1];
			basis->me[5][0] = point[0]*point[1];

			if ( compute == 2) {
				// dpdx
				basis->me[1][1] = 1;
				basis->me[3][1] = 2*point[0];
				basis->me[5][1] = point[1];
				// dpdy
				basis->me[2][2] = 1;
				basis->me[4][2] = 2*point[1];
				basis->me[5][2] = point[0];

			}else if ( compute == 3){


				// dpdx
				basis->me[1][1] = 1;
				basis->me[3][1] = 2*point[0];
				basis->me[5][1] = point[1];
				// dpdy
				basis->me[2][2] = 1;
				basis->me[4][2] = 2*point[1];
				basis->me[5][2] = point[0];


				// d2pdx2
				basis->me[3][3] = 2;
				
				// d2pdxdy
				basis->me[5][4] = 1;

				// d2pdydx
				basis->me[5][5] = 1;

				// d2pdy2
				basis->me[4][6] = 2;

			}else{

			}


		}else if (dim == 3){
			// p(x)
			basis->me[0][0] = 1;
			basis->me[1][0] = point[0];
			basis->me[2][0] = point[1];
			basis->me[3][0] = point[2];

			if ( compute == 2){
				// dpdx
				basis->me[0][1] = 0;
				basis->me[1][1] = 1;
				basis->me[2][1] = 0;
				basis->me[3][1] = 0;
				// dpdy
				basis->me[0][2] = 0;
				basis->me[1][2] = 0;
				basis->me[2][2] = 1;
				basis->me[3][2] = 0;
				// dpdz
				basis->me[0][3] = 0;
				basis->me[1][3] = 0;
				basis->me[2][3] = 0;
				basis->me[3][3] = 1;
			}else {

			}
		}


		}else{
			fprintf(stderr,"dimension not set\n");
		}


		return 0;
}
