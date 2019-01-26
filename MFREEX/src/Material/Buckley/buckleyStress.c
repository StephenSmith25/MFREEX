
#include "Material/Buckley/buckleyStress.h"
#include "determinant.h"

int buckleyStress(state_Buckley * stateNew, state_Buckley * stateOld, VEC * matParams, 
	VEC * critLambdaParams, double dt, int i, int IS_AXI)
{


			int dim = stateNew->F->m;
			double Kb = matParams->ve[9];


			
			// // update to D_n_h and W_n_h
			m_sub(stateNew->F, stateOld->F, stateNew->GRAD_U);
			m_mlt(stateNew->GRAD_U,stateOld->invF,stateNew->H);
			m_add(stateNew->F,stateOld->F,stateNew->delta_F);
			sm_mlt(0.5,stateNew->delta_F,stateNew->delta_F);
			m_inverse_small(stateNew->delta_F,stateNew->invF);
			velocity_grad(stateNew->H, stateNew->L,stateNew->D, stateNew->W,dt,0.500);


			// inverse deformation gradient			
			m_inverse_small(stateNew->F, stateNew->invF);
	




			// Find deviatoric part of the deformation 
			m_copy(stateNew->D,stateNew->Dbar);

			if ( dim == 2)
			{
				stateNew->div_v = stateNew->L->me[0][0] + stateNew->L->me[1][1];

			}else if ( dim == 3)
			{
				stateNew->div_v = stateNew->L->me[0][0] + stateNew->L->me[1][1] + stateNew->L->me[2][2];

			}

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

			//------------------------//
			//      Velocity grad     //
			//------------------------//
			// Find velocity gradient 
			// velocity_grad(stateNew->h, stateNew->L,stateNew->D, stateNew->W,dt,0.500);

			// // div_v
			// if ( dim == 2)
			// {
			// 	stateNew->div_v = stateNew->L->me[0][0] + stateNew->L->me[1][1];

			// }else if ( dim == 3)
			// {
			// 	stateNew->div_v = stateNew->L->me[0][0] + stateNew->L->me[1][1] + stateNew->L->me[2][2];

			// }

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
			//stateNew->Jacobian = stateOld->Jacobian + stateOld->Jacobian*stateNew->div_v*dt;

			stateNew->Jacobian = determinant(stateNew->F);
			/* ------------------------------------------*/
			/* ---------Isochoric   Deformation---------*/
			/* ------------------------------------------*/

			/* Distortional deformation gradient */
			__smlt__(stateNew->F->base, pow(stateNew->Jacobian, -1.00/3.00), stateNew->Fbar->base, dim*dim);
		

			// FIND THE ROTATION TENSOR
			m_mlt(stateOld->V,stateNew->D,stateNew->temp);
			m_mlt(stateNew->D,stateOld->V,stateNew->temp1);
			m_sub(stateNew->temp1,stateNew->temp,stateNew->temp);

			MAT * Vdot = stateNew->Vdot;
			VEC * h = stateNew->h;
			VEC * w = stateNew->w;
			VEC * omega = stateNew->omega;
			VEC * z = stateNew->z;

			// ZERO INPUTS
			v_zero(h);
			v_zero(w);
			v_zero(omega);
			v_zero(z);
			m_zero(stateNew->Omega);


			// zi = eigejk Vjm Dbamk
			for ( int i = 0  ; i < dim ; i++)
			{
				for ( int j = 0 ; j < dim ; j++)
				{
					for ( int k = 0 ; k < dim ; k++)
					{
						int i1 = i+1;
						int j1 = j+1;
						int k1 = k+1;
						double e = ((i1-j1)*(j1-k1)*(k1-i1))/2;
						z->ve[j] += (double)(0.5*e)*stateNew->temp->me[i][k];
						w->ve[j] += (double)(0.5*e)*stateNew->W->me[i][k];

					}
				}

			}


			m_ident(stateNew->temp);

			double traceV = stateOld->V->me[0][0] + stateOld->V->me[1][1] + stateOld->V->me[2][2];
			sm_mlt(traceV,stateNew->temp,stateNew->temp);
			m_sub(stateNew->temp,stateOld->V, stateNew->temp);

			CHfactor(stateNew->temp);
			CHsolve(stateNew->temp, z, h);

			v_add(w,h,omega);




			for ( int i = 0  ; i < dim ; i++)
			{
				for ( int j = 0 ; j < dim ; j++)
				{
					for ( int k = 0 ; k < dim ; k++)
					{
						int i1 = i+1;
						int j1 = j+1;
						int k1 = k+1;
						int e = ((i1-j1)*(j1-k1)*(k1-i1))/2;
						stateNew->Omega->me[i][k] += e*omega->ve[j];

					}
				}

			}


			m_mlt(stateNew->L,stateOld->V,Vdot);
			m_mlt(stateOld->V,stateNew->Omega,stateNew->temp);
			m_sub(Vdot,stateNew->temp,Vdot);

			symmeig(Vdot, stateNew->temp, stateNew->lambdaDot);

			sm_mlt(dt,Vdot,stateNew->V);
			m_add(stateNew->V,stateOld->V,stateNew->V);

			// //Find rotation tensor using Hughes, Winget algorithm
			m_ident(stateNew->temp);
			sm_mlt(0.5*dt,stateNew->Omega,stateNew->temp1);
			m_sub(stateNew->temp,stateNew->temp1,stateNew->temp1);
			m_inverse_small(stateNew->temp1,stateNew->temp);

			m_ident(stateOld->temp);
			sm_mlt(0.5*dt,stateNew->Omega,stateNew->temp1);
			m_add(stateOld->temp,stateNew->temp1,stateNew->temp1);

			m_mlt(stateNew->temp1,stateOld->R,stateOld->temp);
			m_mlt(stateNew->temp,stateOld->temp,stateNew->R);

			// Find co-rotated Dbar and Vdot



	


			// Find delta V
			m_inverse_small(stateOld->V, stateOld->temp);
			m_mlt(stateNew->V,stateOld->temp,stateNew->delta_V);

			m_ident(stateNew->temp);
			m_sub(stateNew->delta_V,stateNew->temp,stateNew->temp);


			sm_mlt(1.00/dt,stateNew->temp,stateNew->temp1);


			// Rotate d to rotation-neutralised configuration 

			mtrm_mlt(stateNew->R,stateNew->temp1,stateNew->temp);
			m_mlt(stateNew->temp,stateNew->R,stateNew->d);


			// v_foutput(stdout,w);
			// v_foutput(stdout,h);
			// v_foutput(stdout,z);
			// v_foutput(stdout,omega);

			// printf("Omega = ");
			// m_foutput(stdout,stateNew->Omega);

			// printf("F = ");
			// m_foutput(stdout,stateNew->F);
			// // printf("R = ");
			// // m_foutput(stdout,stateNew->R);
			// printf("Vdot = ");
			// m_foutput(stdout,Vdot);

			// printf("V = ");
			// m_foutput(stdout,stateNew->V);
			// m_foutput(stdout, stateNew->F);
			//printf("Vdot = ");
			//m_foutput(stdout,Vdot);
		






			// // Polar Decomposition
			// poldec(stateNew->Fbar, stateNew->R, stateNew->U, stateNew->V);

			// Get priciple nominal strain rates 
			tracecatch(
			symmeig(Vdot, stateNew->temp, stateNew->lambdaDot);,
			"Eigen values of V in internalForce");

			tracecatch(
			symmeig(stateNew->Dbar,stateNew->eigVecDBar,stateNew->eigValDBar);,
			"Eigen values of D in internalForce");
			
			PERM * order = px_get(dim);
			stateNew->eigValDBar = v_sort(stateNew->eigValDBar, order);
			px_free(order);




			// // update critical network stretch 
			stateNew->critLambdaBar =lambdaCrit(stateOld->critLambdaBar,stateNew,
			 critLambdaParams, stateNew->temperature, dt, IS_AXI);

			/* ------------------------------------------*/
			/* ------------- --Update Stress--------------*/
			/* ------------------------------------------*/
			/*  Obtain stresses using explicit integration of stress rate */
			buckleyBond(stateNew,stateOld,matParams,dt);
			// Conformational stress
			buckleyConf(stateNew,stateOld, matParams,dt, IS_AXI);

			m_mlt(stateNew->R,stateNew->Sb_R,stateNew->temp);
			mmtr_mlt(stateNew->temp,stateNew->R,stateNew->Sb);


			// m_mlt(stateNew->R,stateNew->Sc_R,stateNew->temp);
			// m_mlt(stateNew->temp,stateNew->R,stateNew->Sc);

			// Rotate back



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
			m_copy(stateNew->Sb_R,stateOld->Sb_R);			
			m_copy(stateNew->Sc_R,stateOld->Sc_R);			
			m_copy(stateNew->sigma,stateOld->sigma);
			stateOld->Jacobian = stateNew->Jacobian;
			stateOld->mSigma = stateNew->mSigma;		
			stateOld->div_v = stateNew->div_v;	
			stateOld->critLambdaBar = stateNew->critLambdaBar;			
			m_copy(stateNew->V,stateOld->V);
			m_copy(stateNew->R,stateOld->R);






	return 0;
}