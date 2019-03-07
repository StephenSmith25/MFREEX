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



MAT * BMAT(MAT * Bmat, shape_function * basis_functions, int dim, int is_axi, double r )
{

	int num_points = basis_functions->dphi->m;


	if (Bmat == MNULL)
	{
		if ( is_axi == 1)
		{
		Bmat = m_get(5,2*num_points);
		}else{
		Bmat = m_get(dim*dim,dim*num_points);	
		}
	}

	if ( is_axi == 1)
	{
		Bmat = m_resize(Bmat,5,2*num_points);
	}else{
		Bmat = m_resize(Bmat,dim*dim,dim*num_points);	
	}

	m_zero(Bmat);

	for ( int i = 0 ; i <num_points ; i++ )
	{

		if ( dim == 2){
		
			Bmat->me[0][2*i] += basis_functions->dphi->me[i][0];
			Bmat->me[1][2*i+1] += basis_functions->dphi->me[i][1];
			Bmat->me[2][2*i] += basis_functions->dphi->me[i][1];
			Bmat->me[3][2*i+1] += basis_functions->dphi->me[i][0];

			if ( is_axi == 1)
			{
				if ( r < 1e-3 )
				{
					Bmat->me[4][2*i] += basis_functions->dphi->me[i][0];
				}else{
					Bmat->me[4][2*i] += basis_functions->phi->ve[i]/r;

				}
				//Bmat->me[4][2*i] += phi_der->me[i][2];
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