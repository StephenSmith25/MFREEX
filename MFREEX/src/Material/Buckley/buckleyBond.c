/*
 * =====================================================================================
 *
 *       Filename:  buckleyBond.c
 *
 *    Description:  Finds the bond stress component of the Buckley model
 *
 *        Version:  1.0
 *        Created:  04/05/18 11:09:36
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  STEPHEN SMITH
 *   Organization:  QUEEN'S UNIVERSITY BELFAST 
 *
 * =====================================================================================
 */
#include "Material/Buckley/buckleyBond.h"
#include <math.h>
int buckleyBond(state_variables * stateNew, state_variables * stateOld , VEC * para, double dt)
{


	// material constants
	double Gb = para->ve[10];
	double H0 = para->ve[7];
	double R = para->ve[8];
	double star_T = para->ve[4];
	double temperature = stateNew->temperature;
	double Vs = para->ve[0];
	double Vp = para->ve[1];
	double star_mu0 = para->ve[2];
	double vogel_T = para->ve[3];
	double Cv = para->ve[6];



	MAT * Sb_n = stateOld->Sb_R;
	MAT * Sb_n_1 = stateNew->Sb_R;
	MAT * d = stateNew->dbar;
	MAT * deltaSb = stateOld->m_temp3;

	m_zero(deltaSb);


	double tauOCT = sqrt (  (1.000/3.00) * contraction(Sb_n,Sb_n)  ) ;
	double alpha_sig = 1;
	
	double sigma_m = stateOld->mSigma;


	// if (sigma_m < 0)
	//  {
	//  	sigma_m = 0;
 // 	}

	if ( tauOCT == 0){
		alpha_sig = 1;
	}else{
		alpha_sig = ( Vs*tauOCT/(2*R*temperature) ) * exp ( -Vp * sigma_m/(R*temperature) ) / ( sinh ( Vs * tauOCT/(2*R*temperature)) ) ; 
	}



	// alpha_s
	double alpha_s = exp ( Cv / ( temperature - vogel_T) - Cv / ( star_T - vogel_T) );	
	// alpha_T
	double alpha_T = exp ( (H0/R) * ( 1/temperature - 1/star_T) );
	

	// tau = tau_s * alpha_s * alpha_T * alpha_sig ;
	double tau = (star_mu0/(2*Gb))*alpha_sig*alpha_s*alpha_T;

	stateNew->tau = tau;

	/* Rotation neutralised method */

	double sbFactor = 1 - exp(-dt/tau);

	// mat1 = 2Gb*D*tau
	sm_mlt(2*Gb*tau,d,stateNew->m_temp1);

	// 2Gb*D*tau - Sb_n
	m_sub(stateNew->m_temp1,Sb_n,stateNew->m_temp2);

	// deltaSb = ( 1- exp(-dt/tau)) * ( 2G*D*tau - Sb_n)
	sm_mlt(sbFactor,stateNew->m_temp2,deltaSb);


	m_add(deltaSb,Sb_n,Sb_n_1);



	return 0;
}


// MAT * mat1 = m_get(3,3) ;
// 	MAT * mat2 = m_get(3,3) ;
// 	MAT * deltaSb = m_get(3,3);


// 	double sbFactor = 1 - exp(-dt/tau);
// 	// mat1 = 2Gb*D*tau
// 	sm_mlt(2*Gb*tau,stateNew->Dbar,mat1);

// 	// 2Gb*D*tau - Sb_n
// 	m_sub(mat1,stateOld->Sb,deltaSb);

// 	// deltaSb = ( 1- exp(-dt/tau)) * ( 2G*D*tau - Sb_n)
// 	// corrotational stress increment
// 	sm_mlt(sbFactor,deltaSb, deltaSb);




// 	// Green rate
// 	//W*s
// 	m_mlt(stateNew->Omega,Sb_n,mat1);
// 	// s*W'
// 	m_mlt(Sb_n,stateNew->Omega,mat2);
// 	m_sub(mat1,mat2,mat1);
// 	sm_mlt(dt,mat1,mat1);
// 	m_add(Sb_n,mat1,mat1);
// 	//Sb* = Sb_n + deltaSb + (WSb_n - Sb_n W)dt
// 	m_add(deltaSb,mat1,Sb_n_1);

	

// 	//TRUESDEL RATE
// 	// m_mlt(stateNew->L,Sb_n,mat1);
// 	// mmtr_mlt(Sb_n, stateNew->L, mat2);

// 	// m_add(mat1,mat2,mat1);

// 	// sm_mlt(stateNew->div_v, Sb_n, mat2);

// 	// m_sub(mat1,mat2,Sb_n_1);

// 	// sm_mlt(dt,mat1,Sb_n_1);

// 	// m_add(deltaSb,Sb_n_1,Sb_n_1);
// 	// m_add(Sb_n,Sb_n_1,Sb_n_1);


// 	//Sb_n_1 = Sb_n + deltaSb + (WSb_n - Sb_n W)dt
// 	// m_mlt(stateNew->W,Sb_n_1,mat1);
// 	// m_mlt(Sb_n_1,stateNew->W,mat2);
// 	// m_sub(mat1,mat2,mat1);	
// 	// sm_mlt(0.5*dt,mat1,mat1);

// 	// m_add(Sb_n_1,mat1,Sb_n_1);


// 	M_FREE(deltaSb); 
// 	M_FREE(mat1);
// 	M_FREE(mat2);


