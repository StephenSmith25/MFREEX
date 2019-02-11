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


	if ( stateNew->div_v >= 0 )
	{


		if ( stateOld->lambdaNMax < stateNew->critLambdaBar){
			gamma_n_1 = gammaV(stateNew,stateOld->lambdaNMax,
				stateNew->critLambdaBar,para, deltaT);
			sm_mlt(1.0000/gamma_n_1,stateOld->Sc_R,Ds);

		}else{
			Ds = m_zero(Ds);
		}
	}else{
		Ds = m_zero(Ds);
	}
	double traceDs = (1.00/3.00)*(Ds->me[0][0] + Ds->me[1][1] + Ds->me[2][2]);

	Ds->me[0][0] = Ds->me[0][0] - traceDs;
	Ds->me[1][1] = Ds->me[1][1] - traceDs;
	Ds->me[2][2] = Ds->me[2][2] - traceDs;

	m_zero(Ds);

	// calculate updated network strains 
	MAT * d_total = stateNew->d;

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


	// update maximum network stre



	// double B1, B2, B3, Bmax = 0;
	// B1 = stateNew->Bbar->me[0][0];
	// B2 = stateNew->Bbar->me[1][1];
	// B3 = stateNew->Bbar->me[2][2];


	// Bmax = max(B1,B2);
	// Bmax = max(Bmax,B3);
	// stateNew->lambdaNMax = sqrt(Bmax);

	//stateNew->lambdaNMax = exp(v_max(eigValB,&index));
	double lambda1 = exp(eigValB->ve[0]);
	double lambda2 = exp(eigValB->ve[1]);
	double lambda3 = exp(eigValB->ve[2]);

	stateNew->lambdaNMax = max(lambda1,lambda2);
	stateNew->lambdaNMax = max(stateNew->lambdaNMax,lambda3);

	m_copy(stateNew->ep_n,stateOld->ep_n);
	stateOld->lambdaNMax = stateNew->lambdaNMax;


	// Find rotation neutralised V_n

	// double lambda_1 = eigValB->ve[0] ;
	// double lambda_2 = eigValB->ve[1];
	// double lambda_3 = eigValB->ve[2];


	// double i1 = lambda_1 + lambda_2 + lambda_3;
	// double i2 = lambda_1*lambda_2 + lambda_1*lambda_3 + lambda_2*lambda_3;
	// double i3 = lambda_1*lambda_2*lambda_3;

	// double D = i1*i2 - i3;










	return 0;
}

// int buckleyConf(state_Buckley * stateNew, state_Buckley * stateOld, VEC * para, double deltaT, int IS_AXI)
// {




// 	MAT * intermediate1 = m_get(3,3);
// 	MAT * intermediate2 = m_get(3,3);
// 	MAT * Dn = m_get(3,3);
// 	MAT * Ds = m_get(3,3);


// 	// approach 1 case 1
// 	MAT * BnCr = m_get(3,3);
// 	MAT * BnDot = m_get(3,3);
// 	MAT * deltaB = m_get(3,3);


// 	VEC * Sc_n_1p = v_get(3);
// 	MAT * eigVecB = m_get(3,3);
// 	VEC * eigValB = v_get(3);
// 	MAT * relSpin = m_get(3,3);
// 	MAT * Sc_n_1 = stateNew->Sc;

// 	double Jacobian = stateNew->Jacobian;
// 	double gamma_n_1;
// 	int index = 0; 
// 	++call_count_2;

// 		if ( stateOld->lambdaNMax < stateNew->critLambdaBar){
// 			gamma_n_1 = gammaV(stateNew,stateOld->lambdaNMax,
// 				stateNew->critLambdaBar,para, deltaT, IS_AXI);
// 			sm_mlt(1.0000/gamma_n_1,stateOld->Sc,Ds);
// 				//m_foutput(stdout,stateNew->Sc);
// 		}else{
// 				// Ds = zero;
// 	}


// 	// //APPROACH 1 CASE 1
// 	// network rate of deformation tensor 
// 	m_sub(stateNew->W,stateNew->Omega,relSpin);

// 	m_sub(stateNew->Dbar,Ds,Dn);



// 	//Fixed material frame rate of deformation
// 	//BnCr = Dn_n*Bn_n + Bn_n*Dn_n;
// 	m_add(Dn,relSpin,stateNew->temp);
// 	m_sub(Dn,relSpin,stateNew->temp1);
// 	m_mlt(stateNew->temp,stateOld->Bbar,intermediate1);
// 	m_mlt(stateOld->Bbar,stateNew->temp1,intermediate2);
// 	m_add(intermediate1,intermediate2,BnCr);

// 	//Material time derivative of B  ( Jaumann )
// 	m_mlt(stateNew->Omega,stateOld->Bbar,intermediate1);
// 	m_mlt(stateOld->Bbar,stateNew->Omega,intermediate2);

// 	m_sub(intermediate1,intermediate2,BnDot);
// 	m_add(BnCr,BnDot,BnDot);
// 	sm_mlt(deltaT,BnDot,deltaB);
// 	// Update B
// 	m_add(stateOld->Bbar,deltaB,stateNew->Bbar);



// 	//Material time derivative of B  ( Truesdell )
// 	// m_mlt(stateNew->L,stateOld->Bbar,intermediate1);
// 	// m_mlt(stateOld->Bbar,stateNew->L,intermediate2);
// 	// m_sub(intermediate1,intermediate2,BnDot);
// 	// m_add(BnCr,BnDot,BnDot);
// 	// sm_mlt(deltaT,BnDot,deltaB);





// 	// APPROACH 2 CASE 3-
// 	// m_sub(stateOld->Dbar,Ds,Dn);
// 	// sm_mlt(deltaT,Dn,intermediate1);
// 	// m_ident(intermediate2);
// 	// m_add(intermediate2,intermediate1,intermediate1);
// 	// m_mlt(intermediate1,stateOld->Fn,stateNew->Fn);

// 	// //Spin components
// 	// m_mlt(stateNew->W,stateOld->Fn,intermediate1);
// 	// sm_mlt(deltaT,intermediate1,intermediate1);
// 	// m_add(stateNew->Fn,intermediate1,stateNew->Fn);



// 	// mmtr_mlt(stateNew->Fn, stateNew->Fn, stateNew->Bbar);

// 	// Eigen Value process

// 	tracecatch(symmeig(stateNew->Bbar,eigVecB,eigValB);,	
// 		"Eigen values of Bbar in Buckley conf");
// 	// find edwards vilgis stress
// 	edwardsVilgis(Sc_n_1p,eigValB,para, Jacobian, stateNew->temperature);
// 	// Find /\, such that /\ = Q^T * A * Q;

// 	for ( int i = 0 ; i < Sc_n_1p->max_dim ; i++)
// 	{
// 		Sc_n_1->me[i][i] = Sc_n_1p->ve[i];
// 	}


// 	// rotate stress to global coordinates 
// 	m_mlt(eigVecB,Sc_n_1,intermediate1);
// 	mmtr_mlt(intermediate1,eigVecB,Sc_n_1);


// 	// update maximum network stre

// 	stateNew->lambdaNMax = sqrt(v_max(eigValB,&index));
// 	m_copy(stateNew->Bbar,stateOld->Bbar);
// 	stateOld->lambdaNMax = stateNew->lambdaNMax;

// 	// Approach 2 case 3
// 	m_copy(stateNew->Fn,stateOld->Fn);

// 	M_FREE(intermediate1);
// 	M_FREE(intermediate2);
// 	M_FREE(Dn);
// 	M_FREE(Ds);



// 	//arppoach 1 case 1
// 	M_FREE(BnCr);
// 	M_FREE(BnDot);
// 	M_FREE(deltaB);

// 	m_free(relSpin);


// 	// eigen routines
// 	M_FREE(eigVecB);
// 	V_FREE(eigValB);
// 	V_FREE(Sc_n_1p); 









// 	return 0;
// }
