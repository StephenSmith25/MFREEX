#include "Deformation/velocity_grad.h"


int velocity_grad(state_variables * stateNew, state_variables * stateOld, 
double delta_t,double alpha)
{

	MAT * D = stateNew->D;
	MAT * W = stateNew->W;
	MAT * L = stateNew->L;

	double div_v;

	// Find gradient of incremental displacmeent ( GRAD delta_U)
	MAT * GRAD_U = stateNew->m_temp4;
	m_sub(stateNew->F, stateOld->F, GRAD_U);

	// Find deformation gradient at n+alpha
	MAT * F_n_a = incremental_deformation_gradient(stateNew, stateOld, alpha);
	MAT * invF_n_a = stateOld->m_temp1;
	m_inverse_small(F_n_a, invF_n_a);

	// Find spatial displacment gradient h_n_a =  ( grad_n+alpha delta_u)
	MAT * h = stateOld->m_temp2;
	m_mlt(GRAD_U,invF_n_a,h);


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


	// update div_v ( divergence of velocity gradient)


	if ( dim == 2)
	{
			stateNew->div_v = stateNew->L->me[0][0] + stateNew->L->me[1][1];

			}else if ( dim == 3)
			{
				stateNew->div_v = stateNew->L->me[0][0] + stateNew->L->me[1][1] + stateNew->L->me[2][2];

			}






	return 0;
}