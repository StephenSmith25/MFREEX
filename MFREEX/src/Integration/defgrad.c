

#include "Integration/defgrad.h"


void get_defgrad(MAT * F, MAT * B, VEC * disp, IVEC * neighbours)
{

	m_zero(F);

	for ( int i = 0 ; i < neighbours->max_dim ; i++)
	{
		int indx = neighbours->ive[i];

		if ( B->m == 4)
		{
			F->me[0][0] += B->me[0][2*i]*disp->ve[2*indx];
			F->me[1][1] += B->me[1][2*i+1]*disp->ve[2*indx+1];
			F->me[0][1] += B->me[2][2*i]*disp->ve[2*indx];
			F->me[1][0] += B->me[3][2*i+1]*disp->ve[2*indx+1];
		}

		if ( B->m == 5)
		{
			F->me[0][0] += B->me[0][2*i]*disp->ve[2*indx];
			F->me[1][1] +=  B->me[1][2*i+1]*disp->ve[2*indx+1];
			F->me[0][1] +=  B->me[2][2*i]*disp->ve[2*indx];
			F->me[1][0] += B->me[3][2*i+1]*disp->ve[2*indx+1];
			F->me[2][2] += B->me[4][2*i]*disp->ve[2*indx];

		}


	}


}
