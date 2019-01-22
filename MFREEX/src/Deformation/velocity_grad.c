#include "Deformation/velocity_grad.h"



int velocity_grad(MAT * h, MAT * L, MAT * D, MAT * W, double delta_t, double alpha)
{


	double fact = (1/(2*delta_t));

	int dim = h->m;


	if ( dim == 2)
	{
		D->me[0][0] = fact*(h->me[0][0] + h->me[0][0] + (1.00-2*alpha)* h->me[0][0]*h->me[0][0]);
		D->me[0][1] = fact*(h->me[0][1] + h->me[1][0] + (1.00-2*alpha)* h->me[1][0]*h->me[0][1]);
		D->me[1][0] = fact*(h->me[1][0] + h->me[0][1] + (1.00-2*alpha)* h->me[0][1]*h->me[1][0]);
		D->me[1][1] = fact*(h->me[1][1] + h->me[1][1] + (1.00-2*alpha)* h->me[1][1]*h->me[1][1]);

		// row
		W->me[0][0] = fact*(h->me[0][0] - h->me[0][0]);
		W->me[0][1] = fact*(h->me[0][1] - h->me[1][0]);
		// row2
		W->me[1][0] = fact*(h->me[1][0] - h->me[0][1]);
		W->me[1][1] = fact*(h->me[1][1] - h->me[1][1]);
	}

	if ( dim == 3)
	{
		// Row 1
		D->me[0][0] = fact*(h->me[0][0] + h->me[0][0] + (1.00-2*alpha)* h->me[0][0]*h->me[0][0]);
		D->me[0][1] = fact*(h->me[0][1] + h->me[1][0] + (1.00-2*alpha)* h->me[1][0]*h->me[0][1]);
		D->me[0][2] = fact*(h->me[0][2] + h->me[2][0] + (1.00-2*alpha)* h->me[2][0]*h->me[0][2]);
		// Row 2
		D->me[1][0] = fact*(h->me[1][0] + h->me[0][1] + (1.00-2*alpha)* h->me[0][1]*h->me[1][0]);
		D->me[1][1] = fact*(h->me[1][1] + h->me[1][1] + (1.00-2*alpha)* h->me[1][1]*h->me[1][1]);
		D->me[1][2] = fact*(h->me[1][2] + h->me[2][1] + (1.00-2*alpha)* h->me[2][1]*h->me[1][2]);
		// Row 3
		D->me[2][0] = fact*(h->me[2][0] + h->me[0][2] + (1.00-2*alpha)* h->me[0][2]*h->me[2][0]);
		D->me[2][1] = fact*(h->me[2][1] + h->me[1][2] + (1.00-2*alpha)* h->me[1][2]*h->me[2][1]);
		D->me[2][2] = fact*(h->me[2][2] + h->me[2][2] + (1.00-2*alpha)* h->me[2][2]*h->me[2][2] );


		W->me[0][0] = fact*(h->me[0][0] - h->me[0][0]);
		W->me[0][1] = fact*(h->me[0][1] - h->me[1][0]);
		W->me[0][2] = fact*(h->me[0][2] - h->me[2][0]);

		W->me[1][0] = fact*(h->me[1][0] - h->me[0][1]);
		W->me[1][1] = fact*(h->me[1][1] - h->me[1][1]);
		W->me[1][2] = fact*(h->me[1][2] - h->me[2][1]);

		W->me[2][0] = fact*(h->me[2][0] - h->me[0][2]);
		W->me[2][1] = fact*(h->me[2][1] - h->me[1][2]);
		W->me[2][2] = fact*(h->me[2][2] - h->me[2][2]);
		
	}


	m_add(D,W,L);

	return 0;
}
