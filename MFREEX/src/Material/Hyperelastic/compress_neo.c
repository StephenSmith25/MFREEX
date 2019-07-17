#include "Material/Hyperelastic/compress_neo.h"
#include "matrix2.h"
#include "matrix.h"

int compress_neo(VEC * stressVoigt, state_variables * state, VEC * params){




	double mu = params->ve[0];
	double lambda = params->ve[1];


	m_zero(state->m_temp1);
	m_zero(state->m_temp2);
	m_zero(state->m_temp3);
	m_zero(state->m_temp4);


	/*  Find C */
	mtrm_mlt(state->F,state->F,state->C);

	// inverse C = m_temp_1
	m_inverse_small(state->C, state->m_temp1);


	// Find indeity 
	m_ident(state->m_temp2);

	// I - c^-1
	m_sub(state->m_temp2,state->m_temp1,state->m_temp2);


	// Find second piola kirchoff 
	sm_mlt(mu,state->m_temp2, state->m_temp2);


	double lnj = log(state->Jacobian);

	sm_mlt(lnj*lambda, state->m_temp1, state->m_temp1);

	m_add(state->m_temp1,state->m_temp2,state->m_temp4);

	m_mlt(state->F,state->m_temp4,state->m_temp3);





	if ( stressVoigt->max_dim == 4)
	{
	stressVoigt->ve[0] = state->m_temp3->me[0][0];
	stressVoigt->ve[1] = state->m_temp3->me[1][1];
	stressVoigt->ve[2] = state->m_temp3->me[0][1];
	stressVoigt->ve[3] = state->m_temp3->me[1][0];
	}else if ( stressVoigt->max_dim == 5)
	{
	stressVoigt->ve[0] = state->m_temp3->me[0][0];
	stressVoigt->ve[1] = state->m_temp3->me[1][1];
	stressVoigt->ve[2] = state->m_temp3->me[0][1];
	stressVoigt->ve[3] = state->m_temp3->me[1][0];	
	stressVoigt->ve[4] = state->m_temp3->me[2][2];	

	}


	// m_temp3 contaisn the first piola kirchoff stress
	sm_mlt(1.000/state->Jacobian,state->m_temp3,state->m_temp3);
	mmtr_mlt(state->m_temp3, state->F, state->sigma);	



return 0;
}
