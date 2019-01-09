
#include "Material/Buckley/buckleyStress.h"



int buckleyStress(state_Buckley * stateNew, state_Buckley * stateOld, VEC * matParams, VEC * critLambdaParams, double dt)
{


			int dim = stateNew->F->m;
			double Kb = matParams->ve[9];


			/* Find Fdot */
			m_mlt(stateNew->F, stateOld->invF, stateNew->delta_F);
			m_sub(stateNew->F, stateOld->F, stateNew->delta_F);

			sm_mlt(1.000/dt,stateNew->delta_F,stateNew->Fdot);
			

			// inverse deformation gradient
			m_inverse(stateNew->F,stateNew->invF);
			poldec(stateNew->delta_F, stateNew->delta_R, stateNew->delta_U, stateNew->delta_V);


			//------------------------//
			//      Velocity grad     //
			//------------------------//
			// Find velocity gradient 
			velocity_grad(stateNew->L, stateNew->D, stateNew->W,stateNew->Fdot,stateNew->invF);

			// div_v
			if ( dim == 2)
			{
				stateNew->div_v = stateNew->L->me[0][0] + stateNew->L->me[1][1];

			}else if ( dim == 3)
			{
				stateNew->div_v = stateNew->L->me[0][0] + stateNew->L->me[1][1] + stateNew->L->me[2][2];

			}

			//------------------------//
			//   Strain increment     //
			//------------------------//


			// sm_mlt(dt, stateOld->D, delta_ep_true);
			// // update total strain
			// m_mlt(R, stateOld->ep_true, temp);
			// mmtr_mlt(temp, R, stateOld->ep_true);
			// m_add(stateOld->ep_true, delta_ep_true, stateNew->ep_true);


			// double delta_ep_vol = 0;

			// if ( dim == 2)
			// {
			// 	delta_ep_vol = (1.00/3.00)*(delta_ep_true->me[0][0] +delta_ep_true->me[1][1]);
			// 	delta_ep_true_vol->me[0][0] = delta_ep_vol;
			// 	delta_ep_true_vol->me[1][1] = delta_ep_vol;

			// }else if ( dim == 3)
			// {
			// 	delta_ep_vol = (1.00/3.00)*(delta_ep_true->me[0][0] +delta_ep_true->me[1][1] +delta_ep_true->me[2][2] );
			// 	delta_ep_true_vol->me[0][0] = delta_ep_vol;
			// 	delta_ep_true_vol->me[1][1] = delta_ep_vol;
			// 	delta_ep_true_vol->me[2][2] = delta_ep_vol;
			// }

			// m_sub(delta_ep_true,delta_ep_true_vol,delta_ep_true_iso);


			// Update Jacobian
			stateNew->Jacobian = stateOld->Jacobian + stateOld->Jacobian*stateNew->div_v*dt;


			/* ------------------------------------------*/
			/* ---------Isochoric   Deformation---------*/
			/* ------------------------------------------*/

			/* Distortional deformation gradient */
			__smlt__(stateNew->F->base, pow(stateNew->Jacobian, -1.00/3.00), stateNew->Fbar->base, dim*dim);
			
			m_copy(stateNew->D,stateNew->Dbar);
			if ( dim == 2)
			{
				stateNew->Dbar->me[0][0] -= (1.00/3.00)*stateNew->div_v;
				stateNew->Dbar->me[1][1] -= (1.00/3.00)*stateNew->div_v;

			}else {
				stateNew->Dbar->me[0][0] -= (1.00/3.00)*stateNew->div_v;
				stateNew->Dbar->me[1][1] -= (1.00/3.00)*stateNew->div_v;
				stateNew->Dbar->me[2][2] -= (1.00/3.00)*stateNew->div_v;

			}

			// Polar Decomposition
			poldec(stateNew->Fbar, stateNew->R, stateNew->U, stateNew->V);

			// Get eigen values of polar decomposition
			tracecatch(
			symmeig(stateNew->V,stateNew->eigVecVBar,stateNew->eigValVBar);,
			"Eigen values of V in internalForce");

			tracecatch(
			symmeig(stateNew->Dbar,stateNew->eigVecDBar,stateNew->eigValDBar);,
			"Eigen values of D in internalForce");
			

			// update critical network stretch 
			stateNew->critLambdaBar =lambdaCrit(stateOld->critLambdaBar,stateNew->eigValVBar, 
				stateNew->eigValDBar, critLambdaParams, stateNew->temperature);


			/* ------------------------------------------*/
			/* ------------- --Update Stress--------------*/
			/* ------------------------------------------*/
			/*  Obtain stresses using explicit integration of stress rate */
			buckleyBond(stateNew,stateOld,matParams,dt);
			// Conformational stress
			buckleyConf(stateNew,stateOld, matParams,dt);
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



			m_copy(stateNew->F,stateOld->F);
			m_copy(stateNew->invF,stateOld->invF);			
			m_copy(stateNew->Fbar,stateOld->Fbar);			
			m_copy(stateNew->Sb,stateOld->Sb);			
			m_copy(stateNew->Sc,stateOld->Sc);			
			m_copy(stateNew->sigma,stateOld->sigma);
			stateOld->mSigma = stateNew->mSigma;			
			stateOld->critLambdaBar = stateNew->critLambdaBar;			




	return 0;
}