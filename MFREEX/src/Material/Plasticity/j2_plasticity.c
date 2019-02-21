/* takes a finite strain increment \delta epsilon = \int R' D R ~= R'D'R dt
 and returns the stress, along with updating the plastic strain parameter gamma 
*/



#include "Material/Plasticity/j2_plasticity.h"


#define ROOT_2_OVER_3 0.81649658092

int j2_plasticity(state_variables * stateNew, state_variables * stateOld, VEC * params, double dt)
{




	// **------------------------- PLASTICITY ROUTINE-------------------------------- * //



	double lambda = params->ve[0];
	double mu = params->ve[1];
	double sigma_yield = params->ve[2];

	// Initial variables
	MAT * d_n_1 = stateNew->d;
	MAT * sigma_n_1_t = stateNew->m_temp1;
	MAT * sigma_n = stateOld->sigma_R;
	MAT * sigma_n_1 = stateNew->sigma_R;


	// temp storage
	MAT * sigma_dot_t = stateNew->m_temp2;
	MAT * temp_1 = stateOld->m_temp1;
	MAT * temp_2 = stateOld->m_temp2;
	MAT * ident = stateOld->m_temp3;
	MAT * delta_sigma = stateOld->m_temp4;
	MAT * Q = stateNew->m_temp4;
	MAT * S_trial = stateNew->m_temp3;

	m_ident(ident);


	// Find trial stress state
	double traceD = d_n_1->me[0][0] + d_n_1->me[1][1] + d_n_1->me[2][2];
	

	// Sigma_dot
	sm_mlt(lambda*traceD,ident,temp_1);
	sm_mlt(2*mu,d_n_1,temp_2);
	m_add(temp_1,temp_2,sigma_dot_t);
	sm_mlt(dt,sigma_dot_t,delta_sigma);

	// Find trial stress
	m_add(sigma_n,delta_sigma,sigma_n_1_t);



	// Find deviatoric trial stress
	double trace_sigma_trial = sigma_n_1_t->me[0][0]+sigma_n_1_t->me[1][1] + sigma_n_1_t->me[2][2];
	m_ident(ident);
	sm_mlt(trace_sigma_trial*(1.00/3.00),ident,ident);

	m_sub(sigma_n_1_t,ident,S_trial);


	double norm_trial_state = sqrt(contraction(S_trial,S_trial));
	// Yield function, find if material is yielding
	// effective stress
	double effective_stress = sqrt(   (3.00/2.00)*contraction(S_trial,S_trial));
	double f = norm_trial_state - ROOT_2_OVER_3*sigma_yield;




	// Check yield condition
	if ( f < 0 )
	{

		// all deformation is elastic
		m_copy(sigma_n_1_t,sigma_n_1);




	} else if ( f >= 0 )
	{


		// Find normal to the yield surface
		sm_mlt((1.00/norm_trial_state),S_trial,Q);

		// Find plastic consistencey parameter
		double delta_gamma = (1.00/(2*mu))*(f);

		sm_mlt(delta_gamma*2*mu,Q,delta_sigma);
		m_sub(sigma_n_1_t,delta_sigma,sigma_n_1);

		// Find increment in plastic strain

			// Find deviatoric trial stress
		double trace_sigma= sigma_n_1->me[0][0]+sigma_n_1->me[1][1] + sigma_n_1->me[2][2];
		m_ident(ident);
		sm_mlt(trace_sigma*(1.00/3.00),ident,ident);

		m_sub(sigma_n_1,ident,S_trial);


		norm_trial_state = sqrt(contraction(S_trial,S_trial));
		// Yield function, find if material is yielding
		// effective stress
		f = norm_trial_state - ROOT_2_OVER_3*sigma_yield;


	}



	// state n+1 ------>> n 
	m_copy(stateNew->sigma_R,stateOld->sigma_R);

	// rotate unrotated stress back into global n+1 frame
	m_mlt(stateNew->R,stateNew->sigma_R,stateNew->m_temp1);
	mmtr_mlt(stateNew->m_temp1,stateNew->R,stateNew->sigma);



	return 0;
}