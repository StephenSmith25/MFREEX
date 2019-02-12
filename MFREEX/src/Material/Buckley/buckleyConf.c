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
#include "dsyevq3.h"
#include "dsyevv3.h"
#include "dsyevh3.h"

#include "symmeig_small.h"

static int call_count_2 = 0;


int buckleyConf(state_variables * stateNew, state_variables * stateOld, 
	VEC * para, double deltaT)
{


	MAT * Ds = stateNew->m_temp1;
	VEC * Sc_n_1p = stateNew->v_temp1;


	// approach 1 case 1


	MAT * eigVecB = stateNew->eigVecVBar;
	VEC * eigValB = stateNew->eigValVBar;
	MAT * Sc_n_1 = stateNew->Sc_R;

	double Jacobian = stateNew->Jacobian;
	double gamma_n_1;
	int index = 0; 
	++call_count_2;


		if ( stateOld->lambdaNMax < stateNew->critLambdaBar){
			gamma_n_1 = gammaV(stateNew,stateOld->lambdaNMax,
				stateNew->critLambdaBar,para, deltaT);
			sm_mlt(1.0000/gamma_n_1,stateOld->Sc_R,Ds);

		}else{
			Ds = m_zero(Ds);
		}
	


	// calculate updated network strains 
	MAT * d_total = stateNew->dbar;

	MAT * D_n = stateNew->m_temp4;

	m_sub(d_total,Ds,D_n);

	MAT * delta_ep = stateNew->m_temp3;
	sm_mlt(deltaT,D_n,delta_ep);
	m_add(stateOld->ep_n,delta_ep,stateNew->ep_n);



	dsyevq3(stateNew->ep_n->me,eigVecB->me,eigValB->ve);

	edwardsVilgis(Sc_n_1p,eigValB,para, Jacobian, stateNew->temperature);


	// Find /\, such that /\ = Q^T * A * Q;
	for ( int i = 0 ; i < Sc_n_1p->max_dim ; i++)
	{
		Sc_n_1->me[i][i] = Sc_n_1p->ve[i];
	}

	MAT * intermediate1 = stateOld->m_temp1;



	// rotate stress to global coordinates 
	m_mlt(eigVecB,Sc_n_1,intermediate1);
	mmtr_mlt(intermediate1,eigVecB,Sc_n_1);


	// update maximum network stretch
	double lambda1 = exp(stateNew->ep_n->me[0][0]);
	double lambda2 = exp(stateNew->ep_n->me[1][1]);
	double lambda3 = exp(stateNew->ep_n->me[2][2]);

	stateNew->lambdaNMax = max(lambda1,lambda2);
	stateNew->lambdaNMax = max(stateNew->lambdaNMax,lambda3);



	// update old variables n+1 ->>> n 
	m_copy(stateNew->ep_n,stateOld->ep_n);
	stateOld->lambdaNMax = stateNew->lambdaNMax;








	return 0;
}
