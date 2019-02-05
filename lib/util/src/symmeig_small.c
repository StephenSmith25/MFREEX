#include "symmeig_small.h"
#include "math.h"

// Need to add ways to catch improperly sized matricies


double TOLERANCE_IS_NORMAL = 1e-10;
double TOLERANCE_IS_B_ZERO = 1e-10;
double TOLERANCE_IS_COL_ZERO = 1e-10;
double TOL_IS_EQUAL= 1e-7;

#include <math.h>



int symmeig_small(const MAT * A, MAT * Q, VEC * Out)
{

	int dim = A->m;





	if ( dim == 3)
	{
		double A11 = A->me[0][0],A12 = A->me[0][1],A13 = A->me[0][2];
		double A21 = A->me[1][0],A22 = A->me[1][1],A23 = A->me[1][2];
		double A31 = A->me[2][0],A32 = A->me[2][1],A33 = A->me[2][2];


		/* ------------------------------------------*/
		/* --------CHECK IF MATRIX IS DIAGONAL--------*/
		/* ------------------------------------------*/


		double sum = 0;
		int n=0;

		for ( int i = 0 ; i < 3 ; i++)
		{
			for ( int j = 0 ; j < 3 ; j++)
			{
				if ( i != j)
				{
					sum += A->me[i][j]; 

				}
			}
		}

		if ( fabs(sum ) < TOLERANCE_IS_NORMAL)
		{

			Out->ve[0] = A->me[0][0];
			Out->ve[1] = A->me[1][1];
			Out->ve[2] = A->me[2][2];



			Q->me[0][0] = 1;
			Q->me[1][1] = 1;
			Q->me[2][2] = 1;

			return 0;

		}


		/* --------------------------------------------------*/
		/* --------Solve cubic characterstic equation--------*/
		/* --------------------------------------------------*/




		double I1 = A11+A22+A33;
		double I2 = A11*A22 + A22*A33 + A11*A33 - A12*A21 - A23*A32 - A13*A31;
		double I3 = -A13*A22*A31 + A12*A23*A31 + A13*A21*A32 - A11*A23*A32 - A12*A21*A33 + A11*A22*A33;

		double b = I2 - (1.00/3.00)*pow(I1,2);

		double c = (-2.00/(27.00))*I1*I1*I1 + I1*I2/3.000 - I3;
		double xA;	
		double lambda_1, lambda_2,lambda_3;
		if ( fabs(b) < TOLERANCE_IS_B_ZERO)
		{
			xA = -pow(c,(1.00/3.00));

			lambda_1 = xA + I1/3;
			lambda_2 = lambda_1;
			lambda_3 = lambda_1;
		} else{
			double m =  2.00*sqrt(-b/3);
			double n = 3.00*c/(m*b);
			double t = atan2(sqrt(1-n*n), n)/3.00;


			int A = 1;
			lambda_1 = m*cos(t + 2*(A-1)*M_PI/3) + I1/3.00;
			A = 2;
			lambda_2 = m*cos(t + 2*(A-1)*M_PI/3) + I1/3.00;
			A=3;
			lambda_3 = m*cos(t + 2*(A-1)*M_PI/3) + I1/3.00;
		}


		Out->ve[0] = lambda_1;
		Out->ve[1] = lambda_2;
		Out->ve[2] = lambda_3;



		/* --------------------------------------------------*/
		/* ----------------Find EigenVectors-----------------*/
		/* --------------------------------------------------*/

		if ( Q != MNULL)
		{




		// Apply Cayley-Hamilton Theorem (A-lambda_1 I)(A-lambda_2 I)(A_3*I) = 0;
		// see https://en.wikipedia.org/wiki/Eigenvalue_algorithm


		// (A-lambda_i )(A-lambda_2)(A-lambda_3) = 0;

		// Find degree of multiplicativity


		// FIND REPEATED EIGENVALUES
			int l1_eq_l2 = 0;
			int l2_eq_l3 = 0;
			int l1_eq_l3 = 0;

			int num_repeated = 0;

			if (( fabs(lambda_1 - lambda_2)/max(lambda_1,lambda_2)) < TOL_IS_EQUAL)
			{
				printf("lambda 1 and lambda 2 are rpeated\n");
				l1_eq_l2 = 1;
				++num_repeated;
			}


			if (( fabs(lambda_1 - lambda_3)/max(lambda_1,lambda_3)) < TOL_IS_EQUAL)
			{
				printf("lambda 1 and lambda 3  are rpeated\n");
				l1_eq_l3 = 1;
				++num_repeated;

			}

			if (( fabs(lambda_2 - lambda_3)/max(lambda_2,lambda_3)) < TOL_IS_EQUAL)
			{
				l2_eq_l3 = 1;
				++num_repeated;

				printf("lambda 2 and lambda 3 are rpeated\n");
			}





			double A_lambda1[3][3];
			double A_lambda2[3][3];
			double A_lambda3[3][3];



			for ( int i = 0 ; i < 3 ; i++)
			{
				for ( int j = 0 ; j < 3 ;j++)
				{
					if ( i == j)
					{
						A_lambda1[i][j] = A->me[i][j] - lambda_1;
						A_lambda2[i][j] = A->me[i][j] - lambda_2;
						A_lambda3[i][j] = A->me[i][j] - lambda_3;


					}else{
						A_lambda1[i][j] = A->me[i][j] ;
						A_lambda2[i][j] = A->me[i][j] ;
						A_lambda3[i][j] = A->me[i][j] ;
					}
				}
			}


			if ( num_repeated == 0)
			{


				double n1[3][3] ;
				double n2[3][3];
				double n3[3][3];
				memset(n1, 0, 9*sizeof(double));
				memset(n2, 0, 9*sizeof(double));
				memset(n3, 0, 9*sizeof(double));


				for ( int i = 0 ; i < 3 ; i++)
				{
					for ( int j = 0 ; j < 3 ;j++)
					{
						for ( int k = 0 ; k < 3 ; k++)
						{

							n1[i][j] += A_lambda2[i][k] * A_lambda3[k][j];
							n2[i][j] += A_lambda1[i][k] * A_lambda3[k][j];
							n3[i][j] += A_lambda1[i][k] * A_lambda2[k][j];

						}
					}
				}


				int V1_FOUND = 0;
				int V2_FOUND = 0;
				int V3_FOUND = 0;

				for ( int i = 0 ; i < 3 ; i++)
				{
					if ( (V1_FOUND != 1))
					{
						double L = sqrt(pow(n1[0][i],2) +pow(n1[1][i],2) + pow(n1[2][i],2));
						if ( L >= TOLERANCE_IS_COL_ZERO)
						{
							Q->me[0][0] = n1[0][i]/L;
							Q->me[1][0] = n1[1][i]/L;
							Q->me[2][0] = n1[2][i]/L;
							V1_FOUND = 1;
						}

					}
					if ( (V2_FOUND != 1))
					{
						double L = sqrt(pow(n2[0][i],2) +pow(n2[1][i],2) + pow(n2[2][i],2));
						if ( L >= TOLERANCE_IS_COL_ZERO)
						{
							Q->me[0][1] = n2[0][i]/L;
							Q->me[1][1] = n2[1][i]/L;
							Q->me[2][1] = n2[2][i]/L;	
							V2_FOUND = 1;
						}
					}
					if ( (V3_FOUND != 1))
					{
						double L = sqrt(pow(n3[0][i],2) +pow(n3[1][i],2) + pow(n3[2][i],2));
						if ( L >= TOLERANCE_IS_COL_ZERO)
						{
							Q->me[0][2] = n3[0][i]/L;
							Q->me[1][2] = n3[1][i]/L;
							Q->me[2][2] = n3[2][i]/L;
							V3_FOUND = 1;
						}
					}
				}




			}else if ( num_repeated == 1)
			{

				if (l1_eq_l2 ==1)
				{

					double n1[3][3];
					double n2[3][3];
					double n3[3][3];
					memset(n1, 0, 9*sizeof(double));
					memset(n2, 0, 9*sizeof(double));
					memset(n3, 0, 9*sizeof(double));


					for ( int i = 0 ; i < 3 ; i++)
					{
						for ( int j = 0 ; j < 3 ;j++)
						{
							for ( int k = 0 ; k < 3 ; k++)
							{

								n1[i][j] += A_lambda2[i][k] * A_lambda3[k][j];
								n2[i][j] += A_lambda1[i][k] * A_lambda3[k][j];
								n3[i][j] += A_lambda1[i][k] * A_lambda2[k][j];

							}
						}
					}



				}else if ( l1_eq_l3 == 1)
				{

					double n1[3][3];
					double n2[3][3];
					double n3[3][3];
					memset(n1, 0, 9*sizeof(double));
					memset(n2, 0, 9*sizeof(double));
					memset(n3, 0, 9*sizeof(double));

					
					for ( int i = 0 ; i < 3 ; i++)
					{
						for ( int j = 0 ; j < 3 ;j++)
						{
							for ( int k = 0 ; k < 3 ; k++)
							{

								n1[i][j] += A_lambda2[i][k] * A_lambda3[k][j];
								n2[i][j] += A_lambda1[i][k] * A_lambda3[k][j];
								n3[i][j] += A_lambda1[i][k] * A_lambda2[k][j];


							}
						}
					}

				}else if ( l2_eq_l3 == 1)
				{

					double n1[3][3];
					double n2[3][3];
					double n3[3][3];
					memset(n1, 0, 9*sizeof(double));
					memset(n2, 0, 9*sizeof(double));
					memset(n3, 0, 9*sizeof(double));

					
					for ( int i = 0 ; i < 3 ; i++)
					{
						for ( int j = 0 ; j < 3 ;j++)
						{
							for ( int k = 0 ; k < 3 ; k++)
							{

								n1[i][j] += A_lambda2[i][k] * A_lambda3[k][j];
								n2[i][j] += A_lambda1[i][k] * A_lambda3[k][j];
								n3[i][j] += A_lambda1[i][k] * A_lambda2[k][j];

							}
						}
					}// END OF LOOP


					printf("testing 123 \n \n");

					int V1_FOUND = 0;
					int V2_FOUND = 0;
					int V3_FOUND = 0;

					printf("n2 = \n");
					for ( int i = 0 ; i < 3 ; i++)
					{
						for ( int j = 0 ; j < 3 ; j++)
						{
							printf("%10.2E ", n2[i][j]);
						}
						printf("\n");
					}


					for ( int i = 0 ; i < 3 ; i++)
					{
						if ( (V1_FOUND != 1))
						{
							double L = sqrt(pow(n1[0][i],2) +pow(n1[1][i],2) + pow(n1[2][i],2));

							printf("L = %lf \n", L);
							if ( L >= TOLERANCE_IS_COL_ZERO)
							{
								Q->me[0][0] = n1[0][i]/L;
								Q->me[1][0] = n1[1][i]/L;
								Q->me[2][0] = n1[2][i]/L;
								V1_FOUND = 1;
							}

						}
						if ( (V2_FOUND != 1))
						{
							double L = sqrt(pow(n2[0][i],2) +pow(n2[1][i],2) + pow(n2[2][i],2));

							printf("L = %10.2E \n",L);
							if ( L >= TOLERANCE_IS_COL_ZERO)
							{
								Q->me[0][1] = n2[0][i]/L;
								Q->me[1][1] = n2[1][i]/L;
								Q->me[2][1] = n2[2][i]/L;	
								V2_FOUND = 1;
							}
						}
						if ( (V3_FOUND != 1))
						{

							double L = sqrt(pow(n3[0][i],2) +pow(n3[1][i],2) + pow(n3[2][i],2));

							if ( L >= TOLERANCE_IS_COL_ZERO)
							{
								Q->me[0][2] = n3[0][i]/L;
								Q->me[1][2] = n3[1][i]/L;
								Q->me[2][2] = n3[2][i]/L;
								V3_FOUND = 1;
							}
						}
					}


				}





			}
		






		}







	}




	return 0;
}