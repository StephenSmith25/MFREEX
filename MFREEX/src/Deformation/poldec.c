#include "Deformation/poldec.h"





int poldec(MAT *F, MAT *R, MAT *U, MAT * V)
{
	int dim =F->m ;


	// allocate memory for routine
	VEC * eigVals = v_get(dim);
	MAT * Q = m_get(dim,dim);
	MAT * P  = m_get(dim,dim);
	MAT * temp = m_get(dim,dim);


	mtrm_mlt(F, F, U);
	symmeig(U,Q , eigVals);


	for ( int i = 0 ; i < dim ; i++)
	{
		P->me[i][i] = sqrt(eigVals->ve[i]);
	}

	m_mlt(Q,P,temp);
	mmtr_mlt(temp, Q, U);
	temp = m_inverse(U, temp);
	m_mlt(F,temp,R);
	mmtr_mlt(F,R,V );


	
	
	V_FREE(eigVals);
	M_FREE(Q);
	M_FREE(P);
	M_FREE(temp);


	return 0;
}