#include "m_inverse_small.h"


// use this function on 2x2 and 3x3 matricies only 

MAT *  
m_inverse_small(MAT *A, MAT * OUT)
{

	if ( A == MNULL)
	{
		fprintf(stderr, "input matrix is NULL\n" );
		return MNULL;
	}
	int dim = A->m;

	if ( A->m != A->n){
		fprintf(stderr, "cannot find inverse of non-square matrix\n" );
		return MNULL;
	}

	if ( OUT == MNULL)
	{
		MAT * OUT = m_get(dim,dim);

	}

	if ( (OUT->m != dim) || (OUT->n != dim))
	{
		 m_resize(OUT,dim,dim);

	}
	if ( dim == 2)
	{

		double A11=A->me[0][0], A12 = A->me[0][1];
		double A21=A->me[1][0], A22 = A->me[1][1];


		// complementary minors
		double C11, C12, C21, C22;
		C11 = A22;
		C12 = A21;
		C21 = A12;
		C22 = A11;

		// determinant
		double det = determinant(A);

		// inverse = 1/det * ( Cofactor(complementray minors) )^T
		OUT->me[0][0] = (1.00/det)*C11;
		OUT->me[0][1] = -(1.00/det)*C21;
		OUT->me[1][0] = -(1.00/det)*C12;
		OUT->me[1][1] = (1.00/det)*C22;







	}else if ( dim == 3)
	{
		double A11=A->me[0][0], A12 = A->me[0][1] , A13 = A->me[0][2];
		double A21=A->me[1][0], A22 = A->me[1][1],  A23 = A->me[1][2];
		double A31=A->me[2][0], A32 = A->me[2][1],  A33 = A->me[2][2];

		// complementary minors
		double C11, C12, C13;
		double C21, C22, C23;
		double C31, C32, C33;

		C11 = A22*A33 - A23*A32;
		C12 = A21*A33 - A31*A23;
		C13 = A21*A32 - A31*A22;

		C21 = A12*A33 - A32*A13;
		C22 = A11*A33 - A31*A13;
		C23 = A11*A32 - A31*A12;

		C31 = A12*A23 - A22*A13;
		C32 = A11*A23 - A21*A13;
		C33 = A11*A22 - A12*A21;


		// determinant
		double det = determinant(A);

		// Cofactors
		double Co11= C11, Co12= -C12, Co13= C13;
		double Co21= -C21, Co22= C22, Co23 = -C23;
		double Co31 = C31, Co32 = -C32, Co33 = C33;

		OUT->me[0][0] = (1.00/det)*Co11;
		OUT->me[0][1] = (1.00/det)*Co21;
		OUT->me[0][2] = (1.00/det)*Co31;

		OUT->me[1][0] = (1.00/det)*Co12;
		OUT->me[1][1] = (1.00/det)*Co22;	
		OUT->me[1][2] = (1.00/det)*Co32;


		OUT->me[2][0] = (1.00/det)*Co13;
		OUT->me[2][1] = (1.00/det)*Co23;
		OUT->me[2][2] = (1.00/det)*Co33;




	}else{
		fprintf(stderr, "Cannot find inverse of matrix greater than 3x3, using m_inverse instead \n");
	}


	return OUT;
}