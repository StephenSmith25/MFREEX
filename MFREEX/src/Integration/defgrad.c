

#include "Integration/defgrad.h"


void get_defgrad(MAT * f, SCNI * scni, VEC * disp){


	MAT * B = scni->B;
	IVEC * neighbours = scni->sfIndex;
	MAT * F_r = scni->F_r;


	double f11 = 0,f22 = 0,f33=0,f12 = 0, f13 = 0, f21 = 0, f23 = 0, f31 = 0, f32=0;
	f11 = 1;
	f22 = 1;
	f33 = 1;

		// Find incremental deformation gradient (f)
	for ( int i = 0 ; i < neighbours->max_dim ; i++)
	{
		int indx = neighbours->ive[i];

		if ( B->m == 4)
		{
			f11 += B->me[0][2*i]*disp->ve[2*indx];
			f22 += B->me[1][2*i+1]*disp->ve[2*indx+1];
			f12 += B->me[2][2*i]*disp->ve[2*indx];
			f21 += B->me[3][2*i+1]*disp->ve[2*indx+1];
		}

		if ( B->m == 5)
		{
			f->me[0][0] += B->me[0][2*i]*disp->ve[2*indx];
			f->me[1][1] +=  B->me[1][2*i+1]*disp->ve[2*indx+1];
			f->me[0][1] +=  B->me[2][2*i]*disp->ve[2*indx];
			f->me[1][0] += B->me[3][2*i+1]*disp->ve[2*indx+1];
			f->me[2][2] += B->me[4][2*i]*disp->ve[2*indx];

		}


	}

	f->me[0][0] = f11*F_r->me[0][0] + f12*F_r->me[1][0];
	f->me[1][1] = f21*F_r->me[0][1] + f22*F_r->me[1][1];
	f->me[0][1] = f11*F_r->me[0][1] + f12*F_r->me[1][1];
	f->me[1][0] = f21*F_r->me[0][0] + f22*F_r->me[1][0];



}

void get_dot_defgrad(MAT * f, SCNI * scni, VEC * velocity){

	MAT * B = scni->B;
	IVEC * neighbours = scni->sfIndex;
	MAT * F_r = scni->F_r;


	double f11 = 0,f22 = 0,f33=0,f12 = 0, f13 = 0, f21 = 0, f23 = 0, f31 = 0, f32=0;


		// Find incremental deformation gradient (f)
	for ( int i = 0 ; i < neighbours->max_dim ; i++)
	{
		int indx = neighbours->ive[i];

		if ( B->m == 4)
		{
			f11 += B->me[0][2*i]*velocity->ve[2*indx];
			f22 += B->me[1][2*i+1]*velocity->ve[2*indx+1];
			f12 += B->me[2][2*i]*velocity->ve[2*indx];
			f21 += B->me[3][2*i+1]*velocity->ve[2*indx+1];
		}

		if ( B->m == 5)
		{
			f->me[0][0] += B->me[0][2*i]*velocity->ve[2*indx];
			f->me[1][1] +=  B->me[1][2*i+1]*velocity->ve[2*indx+1];
			f->me[0][1] +=  B->me[2][2*i]*velocity->ve[2*indx];
			f->me[1][0] += B->me[3][2*i+1]*velocity->ve[2*indx+1];
			f->me[2][2] += B->me[4][2*i]*velocity->ve[2*indx];

		}


	}

	f->me[0][0] = f11*F_r->me[0][0] + f12*F_r->me[1][0];
	f->me[1][1] = f21*F_r->me[0][1] + f22*F_r->me[1][1];
	f->me[0][1] = f11*F_r->me[0][1] + f12*F_r->me[1][1];
	f->me[1][0] = f21*F_r->me[0][0] + f22*F_r->me[1][0];



}
