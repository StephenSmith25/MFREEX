

#include "matop_3x3.h"


MAT	*m_mlt_3x3(const MAT *A, const MAT *B, MAT *OUT)
{
	/* OUT = A B. Every entry is computed: the previous version skipped
	   OUT[i][k] whenever A[i][k] == 0, which is not a valid shortcut
	   (e.g. I * D lost the off-diagonal terms of D). Use temporaries so
	   OUT may alias A or B. */
	Real	**A_v = A->me, **B_v = B->me;
	Real	tmp[3][3];

	for ( int i=0; i<3; i++ )
		for ( int k=0; k<3; k++ )
			tmp[i][k] = A_v[i][0]*B_v[0][k] + A_v[i][1]*B_v[1][k] + A_v[i][2]*B_v[2][k];

	for ( int i=0; i<3; i++ )
		for ( int k=0; k<3; k++ )
			OUT->me[i][k] = tmp[i][k];

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

