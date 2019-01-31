#include "Deformation/rotate_tensor.h"


int rotate_tensor(MAT * A, MAT * Q, MAT * temp, MAT * OUT)
{

	m_mlt(Q,A,temp);
	mtrm_mlt(temp,Q, OUT);

	return 0;
}

int un_rotate_tensor(MAT * A, MAT * Q, MAT * temp, MAT * OUT){


	mtrm_mlt(Q,A,temp);
	m_mlt(temp,Q, OUT);

	return 0;

}

