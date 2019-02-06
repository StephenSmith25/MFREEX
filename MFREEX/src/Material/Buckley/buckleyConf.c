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




	MAT * intermediate1 = stateNew->m_temp1;
	MAT * intermediate2 = stateNew->m_temp2;
	MAT * intermediate3 = stateNew->m_temp3;


	MAT * Dn = stateOld->Dn;
	MAT * Ds = stateOld->m_temp1;
	MAT * relSpin = stateOld->m_temp2;
	VEC * Sc_n_1p = stateNew->v_temp1;


	// approach 1 case 1


	MAT * eigVecB = stateNew->eigVecVBar;
	VEC * eigValB = stateNew->eigValVBar;
	MAT * Sc_n_1 = stateNew->Sc;

	double Jacobian = stateNew->Jacobian;
	double gamma_n_1;
	int index = 0; 
	++call_count_2;


	if ( stateOld->lambdaNMax < stateNew->critLambdaBar){


		gamma_n_1 = gammaV(stateNew,stateOld->lambdaNMax,
			stateNew->critLambdaBar,para, deltaT);
		sm_mlt(1.0000/gamma_n_1,stateOld->Sc,Ds);




	}else{
		Ds = m_zero(Ds);
	}



	// // //APPROACH 1 CASE 1
	// // network rate of deformation tensor 
	m_sub(stateNew->W,stateNew->Omega,relSpin);
	m_sub(stateNew->Dbar,Ds,Dn);


	//Fixed material frame rate of deformation
	//BnCr = Dn_n*Bn_n + Bn_n*Dn_n;
	m_add(Dn,relSpin,intermediate1);
	m_sub(Dn,relSpin,intermediate2);


	m_mlt(intermediate1,stateOld->Bbar,intermediate3);
	m_mlt(stateOld->Bbar,intermediate2,intermediate1);

	MAT * BnCr = stateOld->m_temp1;
	MAT * BnDot = stateOld->m_temp2;
	MAT * deltaB = stateOld->m_temp3;

	m_add(intermediate1,intermediate3,BnCr);



	//Material time derivative of B  ( Green-Naghdi )
	m_mlt(stateNew->Omega,stateOld->Bbar,intermediate1);
	m_mlt(stateOld->Bbar,stateNew->Omega,intermediate2);

	m_sub(intermediate1,intermediate2,BnDot);



	m_add(BnCr,BnDot,BnDot);
	sm_mlt(deltaT,BnDot,deltaB);
	// Update B
	m_add(stateOld->Bbar,deltaB,stateNew->Bbar);



	//Material time derivative of B  ( Truesdell )
	// m_mlt(stateNew->L,stateOld->Bbar,intermediate1);
	// m_mlt(stateOld->Bbar,stateNew->L,intermediate2);
	// m_sub(intermediate1,intermediate2,BnDot);
	// m_add(BnCr,BnDot,BnDot);
	// sm_mlt(deltaT,BnDot,deltaB);



	//APPROACH 2, with spin
	// m_sub(stateNew->Lbar,stateNew->Omega,stateNew->temp);
	// m_sub(stateNew->temp,Ds,stateNew->temp);
	// sm_mlt(deltaT,stateNew->temp,intermediate1);



	// m_sub(stateNew->Dbar,Ds,Dn);
	// sm_mlt(deltaT,Dn,intermediate1);
	// m_ident(intermediate2);
	// m_add(intermediate2,intermediate1,intermediate1);
	// m_mlt(intermediate1,stateOld->Fn,stateNew->Fn);

	// //Spin components
	// m_mlt(stateNew->Omega,stateOld->Fn,intermediate1);
	// sm_mlt(deltaT,intermediate1,intermediate1);
	// m_add(stateNew->Fn,intermediate1,stateNew->Fn);


	// method 2, no rotation
	// m_sub(stateNew->Dbar,Ds,Dn);
	// sm_mlt(deltaT,Dn,intermediate1);
	// m_ident(intermediate2);
	// m_add(intermediate2,intermediate1,intermediate1);
	// m_mlt(intermediate1,stateOld->Fn,stateNew->Fn);





	//mmtr_mlt(stateNew->Fn, stateNew->Fn, stateNew->Bbar);
	//m_copy(stateNew->Fn,stateOld->Fn);
	
	// Eigen Value process

	// tracecatch(symmeig(stateNew->Bbar,eigVecB,eigValB);,	
	// "Eigen values of Bbar in Buckley conf");
	dsyevh3(stateNew->Bbar->me,eigVecB->me,eigValB->ve);


	// 		printf("for loop print\n");

	// 		for ( int i = 0 ; i < 3 ; i++)
	// 		{
	// 			for ( int j = 0 ; j < 3 ; j++)
	// 			{
	// 				printf("%lf ",eigVecB->me[i][j]);
	// 			}
	// 			printf("\n");
	// 		}
	// // find edwards vilgis stress



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
