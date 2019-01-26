#include "Material/Buckley/new_Buckley_State.h"

state_Buckley ** new_Buckley_State(int num_Points, double * temperatures, int is_AXI, int dim)
{

	state_Buckley ** state = malloc(num_Points * sizeof(state_Buckley*));

	int dim_s = dim;
	if ( is_AXI == 1)
	{
		dim_s = 3;
	}

	

	for ( int i = 0 ; i < num_Points; i++){

		state[i] =  malloc(sizeof(state_Buckley));

		state[i]->F = m_get(dim_s,dim_s);
		m_ident(state[i]->F);
		state[i]->Fdot = m_get(dim_s,dim_s);
		state[i]->invF = m_get(dim_s,dim_s);
		m_ident(state[i]->invF);

		state[i]->Fn = m_get(dim_s,dim_s);
		m_ident(state[i]->Fn);

		state[i]->Vdot = m_get(dim_s,dim_s);

		state[i]->delta_F = m_get(dim_s,dim_s);
		state[i]->delta_R = m_get(dim_s,dim_s);
		state[i]->delta_U = m_get(dim_s,dim_s);
		state[i]->delta_V = m_get(dim_s,dim_s);

		state[i]->Fbar = m_get(dim_s,dim_s);
		m_ident(state[i]->Fbar);
		state[i]->Fbardot = m_get(dim_s,dim_s);
		state[i]->invFbar = m_get(dim_s,dim_s);

		state[i]->R = m_get(dim_s,dim_s);
		state[i]->U = m_get(dim_s,dim_s);
		state[i]->V = m_get(dim_s,dim_s);
		m_ident(state[i]->V);
		m_ident(state[i]->R);
	
		state[i]->L = m_get(dim_s,dim_s);
		state[i]->D = m_get(dim_s,dim_s);
		state[i]->W = m_get(dim_s,dim_s);
		state[i]->d = m_get(dim_s,dim_s);

		state[i]->Lbar = m_get(dim_s,dim_s);
		state[i]->Dbar = m_get(dim_s,dim_s);
		state[i]->Bbar = m_get(dim_s,dim_s);
		m_ident(state[i]->Bbar);
		state[i]->Omega = m_get(dim_s,dim_s);
		state[i]->Wbar = m_get(dim_s,dim_s);


		state[i]->eigVecDBar = m_get(dim_s,dim_s);
		state[i]->eigValDBar = v_get(dim_s);
		state[i]->eigVecVBar = m_get(dim_s,dim_s);
		state[i]->eigValVBar = v_get(dim_s);
		state[i]->lambdaDot = v_get(dim_s);


		state[i]->true_strain = m_get(dim_s,dim_s);


		state[i]->Sb = m_get(dim_s,dim_s);
		state[i]->Sc = m_get(dim_s,dim_s);
		state[i]->Sb_R = m_get(dim_s,dim_s);
		state[i]->Sc_R = m_get(dim_s,dim_s);

		state[i]->sigma = m_get(dim_s,dim_s);


		state[i]->lambdaNBar = v_get(dim_s);
		state[i]->critLambdaBar = 100.00;
		state[i]->mSigma = 0;
		state[i]->temperature = temperatures[i]+273.15;

		state[i]->div_v = 0;
		state[i]->lambdaNMax = 1;
		state[i]->Jacobian = 1;

		state[i]->H = m_get(dim_s,dim_s);
		state[i]->GRAD_U = m_get(dim_s,dim_s);

		state[i]->temp = m_get(dim_s,dim_s);
		state[i]->temp1 = m_get(dim_s,dim_s);


		state[i]->w = v_get(dim_s);
		state[i]->h = v_get(dim_s);
		state[i]->z = v_get(dim_s);
		state[i]->omega = v_get(dim_s);


		state[i]->gamma = 0;



	}




	return state;
}
