#include "Bmat/generate_Bmat.h"


MAT * generate_Bmat(MAT * phi_der, int dim, int is_axi, double r )
{

	MAT * Bmat;



	int num_points = phi_der->m;
	// additional row on Bmat
	if ( is_axi == 1)
	{
		Bmat = m_get(5,2*num_points);
	}else{
		Bmat = m_get(dim*dim,dim*num_points);	
	}



	for ( int i = 0 ; i <num_points ; i++ )
	{

		if ( dim == 2){
		
			Bmat->me[0][2*i] += phi_der->me[i][0];
			Bmat->me[1][2*i+1] += phi_der->me[i][1];
			Bmat->me[2][2*i] += phi_der->me[i][1];
			Bmat->me[3][2*i+1] += phi_der->me[i][0];

			if ( is_axi == 1)
			{
				Bmat->me[4][2*i] += phi_der->me[i][2];
			}

		}else if ( dim == 3){

		// not yet implemented 	
		Bmat->me[i][2*i] = 0;
		Bmat->me[i][2*i] = 0;
		Bmat->me[i][2*i] = 0;
		Bmat->me[i][2*i] = 0;

		Bmat->me[i][2*i] = 0;
		Bmat->me[i][2*i] = 0;
		Bmat->me[i][2*i] = 0;
		Bmat->me[i][2*i] = 0;
		Bmat->me[i][2*i] = 0;


		printf("Dimension not yet implemented in Bmat");



		}else{

		}
	}


	return Bmat;

}