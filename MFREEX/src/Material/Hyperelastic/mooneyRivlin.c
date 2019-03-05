
#include "Material/Hyperelastic/mooneyRivlin.h"
#include "matrix2.h"
#include "matrix.h"

int mooneyRivlin(VEC * stressVoigt, state_variables * state, VEC * params){



	/*  Material constants */
	double c1 = params->ve[0];
	double c2 = params->ve[1];
	double lambda = params->ve[2];
	double p0 = c1 + 2*c2;


	m_zero(state->m_temp1);
	m_zero(state->m_temp2);
	m_zero(state->m_temp3);
	m_zero(state->m_temp4);


	/*  Find C */
	mtrm_mlt(state->F,state->F,state->C);


	/*  I1 and I3 */
	double I1 = trace(state->C);
	int dim = state->C->m;
	if ( dim == 2)
	{
		I1 = I1 + 1;
	}

	double I3 = determinant(state->C);


	m_ident(state->m_temp1);
	/*  Stress computation */
	/*  2*c1*d_ij */
	sm_mlt(2*c1,state->m_temp1,state->m_temp2);
	
	/*  I1*ident */
	sm_mlt(I1,state->m_temp1,state->m_temp3);

	// I1 - C
	m_sub(state->m_temp3,state->C,state->m_temp3);
	
	/*  2*c2*(I1 - C) */
	sm_mlt(2*c2,state->m_temp3,state->m_temp3);
		
	/*  Find C^{-1} */
	m_inverse_small(state->C,state->m_temp1);

	/*  Pressure term */
	sm_mlt(-2*(p0-lambda*log(I3)),state->m_temp1,state->m_temp4);

	/*  Find stress */
	m_add(state->m_temp2,state->m_temp3,state->m_temp3);
	m_add(state->m_temp3,state->m_temp4,state->m_temp4);

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





	return 0;
}
