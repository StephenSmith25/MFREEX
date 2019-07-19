#include "Material/Hyperelastic/svk.h"
#include "matrix2.h"
#include "matrix.h"

int svk(VEC * stressVoigt, state_variables * state, VEC * params){

	double lambda = params->ve[0];
	double mu = params->ve[1];



	// Zero temp storage matricies
	m_zero(state->m_temp1);
	m_zero(state->m_temp2);
	m_zero(state->m_temp3);
	m_zero(state->m_temp4);


	/*  Find C */
	mtrm_mlt(state->F,state->F,state->C);


	// Find indeity 
	m_ident(state->m_temp2);


	// Find E
	m_sub(state->C,state->m_temp2,state->m_temp3);
	sm_mlt(0.5000,state->m_temp3,state->m_temp3);

	double traceE = 0;
	for ( int i = 0 ; i < state->m_temp3->m ; i++)
	{
		traceE += state->m_temp3->me[i][i];
	}

	// Find indeity 
	m_ident(state->m_temp2);



	// S = lamba trace E I + 2\mu E
	sm_mlt(traceE*lambda,state->m_temp2,state->m_temp2);
	sm_mlt(2*mu,state->m_temp3,state->m_temp3);
	m_add(state->m_temp3,state->m_temp2,state->m_temp4);


	m_mlt(state->F,state->m_temp4,state->m_temp3);


	// Fill up voigt stress 
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

	// m_temp3 contains the first piola kirchoff stress
	sm_mlt(1.000/state->Jacobian,state->m_temp3,state->m_temp3);
	mmtr_mlt(state->m_temp3, state->F, state->sigma);	



	return 0;
}
