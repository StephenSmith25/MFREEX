

#include "matop_3x3.h"


MAT	*m_mlt_3x3(const MAT *A, const MAT *B, MAT *OUT)
{

	Real	**A_v, **B_v /*, *B_row, *OUT_row, sum, tmp */;

	A_v = A->me;		B_v = B->me;

	m_zero(OUT);
	for ( int i=0; i<3; i++ )
		for ( int k=0; k<3; k++ )
		{
			if ( A_v[i][k] != 0.0 )
			OUT->me[i][k] += A_v[i][0]*B_v[0][k] + A_v[i][1]*B_v[1][k] + A_v[i][2]*B_v[2][k];

		    // if ( A_v[i][k] != 0.0 )
		    //     __mltadd__(OUT->me[i],B_v[k],A_v[i][k],(int)3);
		}


	return OUT;
}

MAT	*m_sub_3x3(const MAT *A, const MAT *B, MAT *OUT)
{


	for ( int i=0; i<3; i++ )
		for ( int k=0; k<3; k++ )
		{
			OUT->me[i][k] = A->me[i][k] - B->me[i][k];
		}

		return OUT;
}

MAT	*m_add_3x3(const MAT *A, const MAT *B, MAT *OUT)
{

	for ( int i=0; i<3; i++ )
		for ( int k=0; k<3; k++ )
		{
			OUT->me[i][k] = A->me[i][k] + B->me[i][k];
		}

		return OUT;
}

