#include "Material/Buckley/new_Buckley_State.h"

state_Buckley ** new_Buckley_State(int num_Points, double * Temperatures, int is_AXI, int dim)
{

	state_Buckley ** state = malloc(num_Points * sizeof(state_Buckley*));

	int dim_s = dim;
	if ( is_AXI = 1)
	{
		dim_s = 3;
	}

	

	for ( int i = 0 ; i < num_Points; i++){
		state[i] =  malloc(sizeof(state_Buckley));
		state[i]->Fbar = m_get(dim_s,dim_s);
		m_ident(state[i]->Fbar); 
		state[i]->F = m_get(dim_s,dim_s);
		m_ident(state[i]->F);

		state[i]->L = m_get(dim_s,dim_s);
		state[i]->D = m_get(dim_s,dim_s);
		state[i]->W = m_get(dim_s,dim_s);

		state[i]->Lbar = m_get(dim_s,dim_s);
		state[i]->Dbar = m_get(dim_s,dim_s);
		state[i]->Bbar = m_get(dim_s,dim_s);
		state[i]->Wbar = m_get(dim_s,dim_s);

		m_ident(state[i]->Bbar);
		state[i]->Wbar = m_get(dim_s,dim_s);
		state[i]->Sc = m_get(dim_s,dim_s);
		state[i]->Sb = m_get(dim_s,dim_s);
		state[i]->eigValDBar = v_get(dim_s);
		state[i]->critLambdaBar = 100.00;
		state[i]->mSigma = 0;
		state[i]->lambdaBar = v_get(dim_s);
		state[i]->lambdaNMax = 1;
		//state[i]->temperature = temperatures[i] ;
		state[i]->Jacobian = 1;

	}




	return state;
}
