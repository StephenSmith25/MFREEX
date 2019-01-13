#include "matmul.h"



int matmul(char TRANS_A, char TRANS_B, int m, int n, double A[],int k, int l, double B[], double C[] )
{

	if (TRANS_A == 'Y')
	{
		if (TRANS_B == 'Y')
		{

		}else{

		}
	}


	if (TRANS_A =='N')
	{
		if (TRANS_B == 'Y')
		{

			if (( m == n) && (k == l)) 
			{	
				if ( m == 2)
				{
					double A11 = 0; double A12 = 0 ; 
					double A21 = 0; double A22 = 0 ;

				}

			}


		}else{
			for ( int i = 0 ; i < m ; i++)
			{
				for ( int j = 0 ; j < l; j++ )
				{

					for ( int r = 0 ; r < k ; r++)
					{

					}
				}





			}
		}
	}


	return 0;
}
