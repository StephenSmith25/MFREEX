
#include "Material/Buckley/buckleyStress.h"
#include "determinant.h"
#include "contraction.h"
#include "Deformation/rotate_tensor.h"
#include "symmeig_small.h"
#include "dsyevc3.h"
#include "dsyevh3.h"
#include "dsyevq3.h"

int 
buckleyStress(state_variables * stateNew, 
	state_variables * stateOld,
	VEC * matParams,
	double dt)

{

			// Dimension of problem
			int dim = stateNew->F->m;

			// get bulk modulus from 
			double Kb = matParams->ve[9];


			/* ------------------------------------------*/
			/* ---------Isochoric   Deformation----------*/
			/* ------------------------------------------*/
			// Find deviatoric part of the deformation 
			m_copy(stateNew->D,stateNew->Dbar);

			if ( dim == 2)
			{
				stateNew->Dbar->me[0][0] += (-1.00/3.00)*stateNew->div_v;
				stateNew->Dbar->me[1][1] += (-1.00/3.00)*stateNew->div_v;

			}else {
				stateNew->Dbar->me[0][0] += (-1.00/3.00)*stateNew->div_v;
				stateNew->Dbar->me[1][1] += (-1.00/3.00)*stateNew->div_v;
				stateNew->Dbar->me[2][2] += (-1.00/3.00)*stateNew->div_v;

			}
			m_add(stateNew->Dbar,stateNew->W,stateNew->Lbar);

			// remove rotation from D to get d
			m_zero(stateNew->m_temp1);
			un_rotate_tensor(stateNew->Dbar, stateNew->R, 
				stateNew->m_temp1, stateNew->dbar);


			// // Find U, and then delta U
			// m_zero(stateNew->m_temp1);
			// mtrm_mlt(stateNew->R, stateNew->V, stateNew->m_temp1);
			// m_mlt(stateNew->m_temp1,stateNew->R,stateNew->U);

			// // Find ubar
			// sm_mlt(pow(stateOld->Jacobian,-1.00/3.00),stateNew->U,stateNew->Ubar);

			// // now find ep_bar
			// m_zero(stateNew->m_temp1);
			// v_zero(stateNew->v_temp1);
			// m_zero(stateNew->m_temp2);
			// m_zero(stateNew->m_temp3);

			// dsyevq3(stateNew->Ubar->me,stateNew->m_temp1->me,stateNew->v_temp1->ve);

			// double ep_1 = log(stateNew->v_temp1->ve[0]);
			// double ep_2 = log(stateNew->v_temp1->ve[1]);
			// double ep_3 = log(stateNew->v_temp1->ve[2]);
			// stateNew->m_temp3->me[0][0] = ep_1;
			// stateNew->m_temp3->me[1][1] = ep_2;
			// stateNew->m_temp3->me[2][2] = ep_3;



			// // rotate strain to global coordinates
			// m_mlt(stateNew->m_temp1,stateNew->m_temp3,stateNew->m_temp2);
			// mmtr_mlt(stateNew->m_temp2,stateNew->m_temp1,stateNew->EP_bar);

			// // Find deltaU
			// // m_zero(stateNew->m_temp2);
			// // m_zero(stateNew->m_temp2);

			// m_sub(stateNew->EP_bar,stateOld->EP_bar,stateNew->delta_ep_bar);


			// m_inverse_small(stateOld->Ubar,stateNew->m_temp2);
			// m_mlt(stateNew->Ubar,stateNew->m_temp2,stateNew->delta_Ubar);

			//m_zero(stateNew->m_temp1);
			// m_mlt(stateOld->R, stateNew->delta_Ubar, stateNew->m_temp1);
			// mmtr_mlt(stateNew->m_temp1,stateOld->R,stateNew->delta_Ubar);


			// // delta ep = ln(delta U ) approx U-I
			// MAT * ident = stateOld->m_temp4;
			// m_ident(ident);
			// m_sub(stateNew->delta_Ubar,ident,stateNew->delta_ep_bar);
	
			// Find Vbar dot
			// m_zero(stateNew->m_temp1);
			// m_zero(stateNew->m_temp2);

			// sm_mlt(pow(stateOld->Jacobian,-1.00/3.00)/dt,stateOld->V,stateNew->m_temp1);
			// sm_mlt(pow(stateNew->Jacobian,-1.00/3.00)/dt,stateNew->V,stateNew->m_temp2);
			// m_sub(stateNew->m_temp2,stateNew->m_temp1,stateNew->Vdot);

			// m_zero(stateNew->m_temp1);
			// mtrm_mlt(stateNew->R,stateNew->Vdot,stateNew->m_temp1);
			// m_mlt(stateNew->m_temp1,stateNew->R,stateNew->Vdot);


			/* ------------------------------------------*/
			/* ----------Eigen value routines------------*/
			/* ------------------------------------------*/

			sm_mlt(pow(stateNew->Jacobian,-1.00/3.00),stateNew->Vdot,stateNew->Vdot);

			dsyevc3(stateNew->Vdot->me, stateNew->lambdaDot->ve);
			// sm_mlt(pow(stateOld->Jacobian,-1.00/3.00)/dt,stateOld->V,stateNew->m_temp1);

			// update critical network stretch 
			stateNew->critLambdaBar =lambdaCrit(stateOld->critLambdaBar,stateNew,
				matParams, stateNew->temperature, dt);

			//stateNew->critLambdaBar = 1.5;


			/* ------------------------------------------*/
			/* ------------- --Update Stress--------------*/
			/* ------------------------------------------*/

			/* Bond stress */
			buckleyBond(stateNew,stateOld,matParams,dt);
			// Conformational stress
			buckleyConf(stateNew,stateOld, matParams,dt);

			// rotate unrotated stress back into n+1 configuration
			m_zero(stateNew->m_temp1);
			m_mlt(stateNew->R,stateNew->Sb_R,stateNew->m_temp1);
			mmtr_mlt(stateNew->m_temp1,stateNew->R,stateNew->Sb);


			// rotate unrotated stress back into n+1 configuration
			m_zero(stateNew->m_temp1);
			m_mlt(stateNew->R,stateNew->Sc_R,stateNew->m_temp1);
			mmtr_mlt(stateNew->m_temp1,stateNew->R,stateNew->Sc);
			//sm_mlt((1.00)/stateNew->Jacobian,stateNew->Sc,stateNew->Sc);

			// Hydrostatic stress
			stateNew->mSigma = log(stateNew->Jacobian)*Kb;

			m_add(stateNew->Sb, stateNew->Sc, stateNew->sigma);

			if ( dim == 2)
			{
				stateNew->sigma->me[0][0] += stateNew->mSigma;
				stateNew->sigma->me[1][1] += stateNew->mSigma;



			}else{
				stateNew->sigma->me[0][0] += stateNew->mSigma;
				stateNew->sigma->me[1][1] += stateNew->mSigma;	
				stateNew->sigma->me[2][2] += stateNew->mSigma;	

			}
		


			// UPDATE PREVIOUS TIME STEP VARAIBLES
			m_copy(stateNew->Sb,stateOld->Sb);		
			m_copy(stateNew->EP_bar,stateOld->EP_bar);	
			m_copy(stateNew->Sc,stateOld->Sc);	
			m_copy(stateNew->Sb_R,stateOld->Sb_R);			
			m_copy(stateNew->Sc_R,stateOld->Sc_R);	
			stateOld->mSigma = stateNew->mSigma;		
			stateOld->critLambdaBar = stateNew->critLambdaBar;			
				




	return 0;
}



			// // Rotate stres
			// mtrm_mlt(stateNew->R,stateNew->sigma,stateNew->temp);
			// m_mlt(stateNew->temp,stateNew->R,stateNew->sig_R);


			// //Find change in stress
			// m_sub(stateNew->sig_R,stateOld->sig_R,stateNew->delta_sig);

		

			// // Find increment in deviatoric strain
			// m_copy(stateNew->d,stateNew->delta_ep_dev);

			// double sigma_vol = stateNew->mSigma - stateOld->mSigma;
			// double ep_vol = (1.00/3.00)*stateOld->Jacobian*stateNew->div_v*dt;


			// // Find stress change
			// m_ident(stateNew->temp);
			// sm_mlt((1.000/3.00)*sigma_vol,stateNew->temp,stateNew->delta_sig_vol);
			// m_sub(stateNew->delta_sig,stateNew->delta_sig_vol,stateNew->delta_sig_dev);


			// double eij_eij = contraction(stateNew->delta_ep_dev, stateNew->delta_ep_dev);
			// double sij_eij = contraction(stateNew->delta_sig_dev, stateNew->delta_ep_dev);
			// double mu;
			// double lambda;
			// double kappa; 


			// double lambda_0 = stateNew->lambda_0;
			// double mu_0 = stateNew->mu_0;



			// if ( fabs(ep_vol) > pow(10,-6))
			// {

			// 	if ( dt*dt*eij_eij > pow(10,-12))
			// 	{
			// 		mu = (1.00/2.00) * ( sij_eij / (dt * eij_eij));
			// 		kappa = 1.8e9;
			// 		lambda = kappa - (2.00/3.00)*mu;


			// 	}else{

			// 		kappa = 1.8e9;
			// 		lambda = lambda_0;
			// 		mu = (1.00/4.00) * ( 3*(lambda_0 + 2*mu_0) - 3*kappa );

			// 	}


			// }else{

			// 	if ( dt*dt*eij_eij > pow(10,-12))
			// 	{
			// 		mu = (1.00/2.00) * ( sij_eij / (dt * eij_eij));
			// 		lambda = lambda_0 + 2*mu_0 - 2*mu;

			// 	}else{
			// 		lambda = lambda_0;
			// 		mu = mu_0;
			// 	}

			// }


			// stateNew->lambda = fabs(lambda);
			// stateNew->mu = fabs(mu);

