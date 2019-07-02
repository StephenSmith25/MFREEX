

#include "Integration/defgrad.h"


void get_defgrad(MAT * f, MAT * B, IVEC * neighbours, MAT * F_r, VEC * disp){


	double f11 = 1,f22 = 1,f33=1, f12 = 0, f13 = 0, f21 = 0, f23 = 0, f31 = 0, f32=0;




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
			f11+= B->me[0][2*i]*disp->ve[2*indx];
			f22+= B->me[1][2*i+1]*disp->ve[2*indx+1];
			f12+= B->me[2][2*i]*disp->ve[2*indx];
			f21+= B->me[3][2*i+1]*disp->ve[2*indx+1];
			f33+= B->me[4][2*i]*disp->ve[2*indx];

		}


	}

	if (B->m == 4){

		f->me[0][0] = f11*F_r->me[0][0] + f12*F_r->me[1][0];
		f->me[1][1] = f21*F_r->me[0][1] + f22*F_r->me[1][1];
		f->me[0][1] = f11*F_r->me[0][1] + f12*F_r->me[1][1];
		f->me[1][0] = f21*F_r->me[0][0] + f22*F_r->me[1][0];
	} else if ( B->m == 5)
	{
	
		f->me[0][0] = f11*F_r->me[0][0] + f12*F_r->me[1][0];
		f->me[1][1] = f21*F_r->me[0][1] + f22*F_r->me[1][1];
		f->me[0][1] = f11*F_r->me[0][1] + f12*F_r->me[1][1];
		f->me[1][0] = f21*F_r->me[0][0] + f22*F_r->me[1][0];
		f->me[2][2] = f33*F_r->me[2][2];
	}



	return;


}

void defgrad(MAT * f, MAT * B, IVEC * neighbours, VEC * disp){



	m_ident(f);


		// Find incremental deformation gradient (f)
	for ( int i = 0 ; i < neighbours->max_dim ; i++)
	{
		int indx = neighbours->ive[i];

		if ( B->m == 4)
		{
			f->me[0][0] += B->me[0][2*i]*disp->ve[2*indx];
			f->me[1][1] += B->me[1][2*i+1]*disp->ve[2*indx+1];
			f->me[0][1] += B->me[2][2*i]*disp->ve[2*indx];
			f->me[1][0] += B->me[3][2*i+1]*disp->ve[2*indx+1];
		}

		if ( B->m == 5)
		{
			f->me[0][0]+= B->me[0][2*i]*disp->ve[2*indx];
			f->me[1][1] += B->me[1][2*i+1]*disp->ve[2*indx+1];
			f->me[0][1]+= B->me[2][2*i]*disp->ve[2*indx];
			f->me[1][0]+= B->me[3][2*i+1]*disp->ve[2*indx+1];
			f->me[2][2]+= B->me[4][2*i]*disp->ve[2*indx];

		}


	}


	return;


}
void defgrad_m(MAT * f, MAT * B, IVEC * neighbours, int num_neighbours, VEC * disp){



	m_ident(f);

	// m_foutput(stdout, B);
	// iv_foutput(stdout,neighbours);

	// Find incremental deformation gradient (f)
	for ( int i = 0 ; i < num_neighbours ; i++)
	{
		int indx = neighbours->ive[i];

		if ( B->m == 4)
		{

			f->me[0][0] += B->me[0][2*i]*disp->ve[2*indx];
			f->me[1][1] += B->me[1][2*i+1]*disp->ve[2*indx+1];
			f->me[0][1] += B->me[2][2*i]*disp->ve[2*indx];
			f->me[1][0] += B->me[3][2*i+1]*disp->ve[2*indx+1];
		}

		if ( B->m == 5)
		{
			f->me[0][0] += B->me[0][2*i]*disp->ve[2*indx];
			f->me[1][1] += B->me[1][2*i+1]*disp->ve[2*indx+1];
			f->me[0][1] += B->me[2][2*i]*disp->ve[2*indx];
			f->me[1][0] += B->me[3][2*i+1]*disp->ve[2*indx+1];
			f->me[2][2] += B->me[4][2*i]*disp->ve[2*indx];

		}


	}

	return;


}

void get_dot_defgrad(MAT * f,MAT * B,IVEC * neighbours, MAT * F_r, VEC * velocity) {


	m_zero(f);
	double f11 = 0,f22 = 0,f33=0,f12 = 0, f13 = 0, f21 = 0, f23 = 0, f31 = 0, f32=0;


	// Find incremental dot deformation gradient (f)
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
			f11 += B->me[0][2*i]*velocity->ve[2*indx];
			f22 += B->me[1][2*i+1]*velocity->ve[2*indx+1];
			f12 += B->me[2][2*i]*velocity->ve[2*indx];
			f21 += B->me[3][2*i+1]*velocity->ve[2*indx+1];
			f33 += B->me[4][2*i]*velocity->ve[2*indx];

		}


	}

	if ( B->m == 4)
	{
		f->me[0][0] = f11*F_r->me[0][0] + f12*F_r->me[1][0];
		f->me[1][1] = f21*F_r->me[0][1] + f22*F_r->me[1][1];
		f->me[0][1] = f11*F_r->me[0][1] + f12*F_r->me[1][1];
		f->me[1][0] = f21*F_r->me[0][0] + f22*F_r->me[1][0];
	}

	if ( B->m == 5)
	{
		f->me[0][0] = f11*F_r->me[0][0] + f12*F_r->me[1][0];
		f->me[1][1] = f21*F_r->me[0][1] + f22*F_r->me[1][1];
		f->me[0][1] = f11*F_r->me[0][1] + f12*F_r->me[1][1];
		f->me[1][0] = f21*F_r->me[0][0] + f22*F_r->me[1][0];
		f->me[2][2] = f33*F_r->me[2][2];

	}
}


void get_grad_u(MAT * f, MAT * B, IVEC * neighbours, MAT * F_r, VEC * disp){


	double f11 = 0,f22 = 0,f33=0, f12 = 0, f13 = 0, f21 = 0, f23 = 0, f31 = 0, f32=0;




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
			f11+= B->me[0][2*i]*disp->ve[2*indx];
			f22+= B->me[1][2*i+1]*disp->ve[2*indx+1];
			f12+= B->me[2][2*i]*disp->ve[2*indx];
			f21+= B->me[3][2*i+1]*disp->ve[2*indx+1];
			f33+= B->me[4][2*i]*disp->ve[2*indx];

		}


	}

	if (B->m == 4){

		f->me[0][0] = f11*F_r->me[0][0] + f12*F_r->me[1][0];
		f->me[1][1] = f21*F_r->me[0][1] + f22*F_r->me[1][1];
		f->me[0][1] = f11*F_r->me[0][1] + f12*F_r->me[1][1];
		f->me[1][0] = f21*F_r->me[0][0] + f22*F_r->me[1][0];
	} else if ( B->m == 5)
	{
	
		f->me[0][0] = f11*F_r->me[0][0] + f12*F_r->me[1][0];
		f->me[1][1] = f21*F_r->me[0][1] + f22*F_r->me[1][1];
		f->me[0][1] = f11*F_r->me[0][1] + f12*F_r->me[1][1];
		f->me[1][0] = f21*F_r->me[0][0] + f22*F_r->me[1][0];
		f->me[2][2] = f33*F_r->me[2][2];
	}



	return;


}