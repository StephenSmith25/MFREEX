/*
 * =====================================================================================
 *
 *       Filename:  buckleyConf.c
 *
 *    Description:  Find the conformational stress in the buckley model
 *
 *        Version:  1.0
 *        Created:  12/05/18 14:21:35
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include "Material/Buckley/buckleyConf.h"




static int call_count_2 = 0;

int buckleyConf(MAT * Sc_n_1, state_Buckley * stateOld, VEC * para,double deltaT,double Jacobian){




	MAT * intermediate1 = m_get(3,3);
	MAT * intermediate2 = m_get(3,3);
	MAT * Dn = m_get(3,3);
	MAT * Ds = m_get(3,3);
	MAT * BnCr = m_get(3,3);
	MAT * BnDot = m_get(3,3);
	MAT * deltaB = m_get(3,3);
	VEC * Sc_n_1p = v_get(3);
	MAT * eigVecB = m_get(3,3);
	VEC * eigValB = v_get(3);

	int index = 0; 
	++call_count_2;

	// conformational stress

	//if ( stateOld->eigValDBar->ve[1] > 0 ){
	if (stateOld->div_v > 0 ){ 

	if ( stateOld->lambdaNMax < stateOld->critLambdaBar){
		double gamma_n_1 = gammaV(stateOld->lambdaBar,stateOld->lambdaNMax,stateOld->critLambdaBar,stateOld->eigValDBar,para);
		sm_mlt(1.0000/gamma_n_1,stateOld->Sc,Ds);

		//m_foutput(stdout,stateOld->Sc);
	}else{
		m_zero(Ds);
	}
	}


	// network rate of deformation tensor 
	m_sub(stateOld->Dbar,Ds,Dn);
	// BnCr = Dn_n*Bn_n + Bn_n*Dn_n; 
	m_mlt(Dn,stateOld->Bbar,intermediate1);
	m_mlt(stateOld->Bbar,Dn,intermediate2);
	m_add(intermediate1,intermediate2,BnCr);
	// Bndot = BnCr + W*B - B*W; 
	m_mlt(stateOld->Wbar,stateOld->Bbar,intermediate1);
	m_mlt(stateOld->Bbar,stateOld->Wbar,intermediate2);
	m_add(BnCr,intermediate1,BnDot);
	m_sub(BnDot,intermediate2,BnDot);
	// deltaB = deltaT * BnDot
	sm_mlt(deltaT,BnDot,deltaB);
	// Bn_n_1 = Bn_n + deltaB;
	m_add(stateOld->Bbar,deltaB,stateOld->Bbar);
	symmeig(stateOld->Bbar,eigVecB,eigValB);	
	// find edwards vilgis stress
	edwardsVilgis(Sc_n_1p,eigValB,para, Jacobian);
	// Find /\, such that /\ = Q^T * A * Q;

	for ( int i = 0 ; i < Sc_n_1p->max_dim ; i++)
	{
		Sc_n_1->me[i][i] = Sc_n_1p->ve[i];
	}


	// rotate stress to global coordinates 
	m_mlt(eigVecB,Sc_n_1,intermediate1);
	mmtr_mlt(intermediate1,eigVecB,Sc_n_1);


	// update maximum network stre

	stateOld->lambdaNMax = sqrt(v_max(eigValB,&index));



	M_FREE(intermediate1);
	M_FREE(intermediate2);
	M_FREE(Dn);
	M_FREE(Ds);
	M_FREE(BnCr);
	M_FREE(BnDot);
	M_FREE(deltaB);
	M_FREE(eigVecB);




	V_FREE(eigValB);
	V_FREE(Sc_n_1p); 









	return 0;
}

