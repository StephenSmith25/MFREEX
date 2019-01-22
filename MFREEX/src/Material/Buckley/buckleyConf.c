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

int buckleyConf(state_Buckley * stateNew, state_Buckley * stateOld, VEC * para, double deltaT)
{




	MAT * intermediate1 = m_get(3,3);
	MAT * intermediate2 = m_get(3,3);
	MAT * Dn = m_get(3,3);
	MAT * Ds = m_get(3,3);


	// approach 1 case 1
	MAT * BnCr = m_get(3,3);
	MAT * BnDot = m_get(3,3);
	MAT * deltaB = m_get(3,3);


	VEC * Sc_n_1p = v_get(3);
	MAT * eigVecB = m_get(3,3);
	VEC * eigValB = v_get(3);

	MAT * Sc_n_1 = stateNew->Sc;

	double Jacobian = stateNew->Jacobian;
	double gamma_n_1;
	int index = 0; 
	++call_count_2;


	if ( stateNew->div_v >= 0){
		if ( stateOld->lambdaNMax < stateNew->critLambdaBar){
			gamma_n_1 = gammaV(stateNew,stateOld->lambdaNMax,
			stateNew->critLambdaBar,para, deltaT);
			sm_mlt(1.0000/gamma_n_1,stateOld->Sc,Ds);
			//m_foutput(stdout,stateNew->Sc);
		}else{
			// Ds = zero;
		}
	}


	

	// //APPROACH 1 CASE 1
	// network rate of deformation tensor 
	m_sub(stateOld->Dbar,Ds,Dn);


	//Fixed material frame rate of deformation
	//BnCr = Dn_n*Bn_n + Bn_n*Dn_n;
	m_mlt(Dn,stateOld->Bbar,intermediate1);
	m_mlt(stateOld->Bbar,Dn,intermediate2);
	m_add(intermediate1,intermediate2,BnCr);

	//Material time derivative of B  ( Jaumann )
	m_mlt(stateOld->W,stateOld->Bbar,intermediate1);
	mmtr_mlt(stateOld->Bbar,stateOld->W,intermediate2);
	m_add(intermediate1,intermediate2,BnDot);
	m_add(BnCr,BnDot,BnDot);
	sm_mlt(deltaT,BnDot,deltaB);
	// Update B
	m_add(stateOld->Bbar,deltaB,stateNew->Bbar);


	//Material time derivative of B  ( Truesdell )
	m_mlt(stateNew->L,stateOld->Bbar,intermediate1);
	m_mlt(stateOld->Bbar,stateNew->L,intermediate2);
	m_sub(intermediate1,intermediate2,BnDot);
	m_add(BnCr,BnDot,BnDot);
	sm_mlt(deltaT,BnDot,deltaB);





	// APPROACH 2 CASE 3-
	// m_sub(stateOld->Dbar,Ds,Dn);
	// sm_mlt(deltaT,Dn,intermediate1);
	// m_ident(intermediate2);
	// m_add(intermediate2,intermediate1,intermediate1);
	// m_mlt(intermediate1,stateOld->Fn,stateNew->Fn);

	// //Spin components
	// m_mlt(stateNew->W,stateOld->Fn,intermediate1);
	// sm_mlt(deltaT,intermediate1,intermediate1);
	// m_add(stateNew->Fn,intermediate1,stateNew->Fn);



	// mmtr_mlt(stateNew->Fn, stateNew->Fn, stateNew->Bbar);


	// Eigen Value process
	symmeig(stateNew->Bbar,eigVecB,eigValB);	
	// find edwards vilgis stress
	edwardsVilgis(Sc_n_1p,eigValB,para, Jacobian, stateNew->temperature);
	// Find /\, such that /\ = Q^T * A * Q;

	for ( int i = 0 ; i < Sc_n_1p->max_dim ; i++)
	{
		Sc_n_1->me[i][i] = Sc_n_1p->ve[i];
	}


	// rotate stress to global coordinates 
	m_mlt(eigVecB,Sc_n_1,intermediate1);
	mmtr_mlt(intermediate1,eigVecB,Sc_n_1);


	// update maximum network stre

	stateNew->lambdaNMax = sqrt(v_max(eigValB,&index));
	m_copy(stateNew->Bbar,stateOld->Bbar);
	stateOld->lambdaNMax = stateNew->lambdaNMax;

	// Approach 2 case 3
	m_copy(stateNew->Fn,stateOld->Fn);

	M_FREE(intermediate1);
	M_FREE(intermediate2);
	M_FREE(Dn);
	M_FREE(Ds);



	//arppoach 1 case 1
	M_FREE(BnCr);
	M_FREE(BnDot);
	M_FREE(deltaB);



	// eigen routines
	M_FREE(eigVecB);
	V_FREE(eigValB);
	V_FREE(Sc_n_1p); 









	return 0;
}

