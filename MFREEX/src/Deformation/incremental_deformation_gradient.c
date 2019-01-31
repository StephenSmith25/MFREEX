#include "Deformation/incremental_deformation_gradient.h"




MAT *  
incremental_deformation_gradient(state_variables * stateNew, state_variables * stateOld,
double alpha)
{


	MAT * F_n_a = stateNew->m_temp3;
	// alpha F_n_1
	sm_mlt(alpha,stateNew->F,stateNew->m_temp1);
	// (1-alpha)F_n
	sm_mlt(1-alpha,stateOld->F,stateNew->m_temp2);
	m_add(stateNew->m_temp1,stateNew->m_temp2,F_n_a);



 	return F_n_a;
 }