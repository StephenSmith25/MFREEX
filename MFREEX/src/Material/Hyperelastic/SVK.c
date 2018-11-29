/*
 * =====================================================================================
 *
 *       Filename:  SVK.c
 *
 *    Description:  Saint Venant-Kirchhoff Material model 
 *
 *        Version:  1.0
 *        Created:  19/02/18 11:45:35
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  STEPHEN SMITH
 *   Organization:  QUEEN'S UNIVERSITY BELFAST
 *
 * =====================================================================================
 */
#include "Material/Hyperelastic/SVK.h"

int SVK(VEC * stressVoigt, MAT * defGrad, VEC * params){



	v_zero(stressVoigt);

	MAT * strainGL = m_get(defGrad->m,defGrad->m);
	MAT * ident = m_get(defGrad->m,defGrad->m);
	MAT * term1 = m_get(defGrad->m,defGrad->m);
	MAT * term2 = m_get(defGrad->m,defGrad->m);
	MAT * stress1PKF = m_get(defGrad->m,defGrad->m);
	MAT * stress2PKF = m_get(defGrad->m,defGrad->m);


	m_ident(ident);


	mtrm_mlt(defGrad,defGrad,strainGL);
	
	m_sub(strainGL,ident,strainGL);
	sm_mlt(0.5000,strainGL,strainGL);
	double traceE = trace(strainGL);




	sm_mlt(params->ve[0]*traceE,ident,term1);
	sm_mlt(2*params->ve[1],strainGL,term2);

	m_add(term1,term2,stress2PKF);

	m_mlt(defGrad,stress2PKF,stress1PKF);



	// if ( AXI == 1){
	// stressVoigt->ve[4] = stress1PKF->me[2][2];
	// }
	stressVoigt->ve[0] = stress1PKF->me[0][0];
	stressVoigt->ve[1] = stress1PKF->me[1][1];
	stressVoigt->ve[2] = stress1PKF->me[0][1];
	stressVoigt->ve[3] = stress1PKF->me[1][0];











	m_free(strainGL);
	m_free(ident);
	m_free(term1);
	m_free(term2);
	m_free(stress1PKF);
	m_free(stress2PKF);



	return 0;







}

