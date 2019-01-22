
#include "Material/Buckley/buckleyStress.h"
#include "determinant.h"

int buckleyStress(state_Buckley * stateNew, state_Buckley * stateOld, VEC * matParams, VEC * critLambdaParams, double dt, int i )
{


			int dim = stateNew->F->m;
			double Kb = matParams->ve[9];


			// /* Find Fdot */
			// m_sub(stateNew->F, stateOld->F, stateNew->delta_F);

			// sm_mlt(1.000/dt,stateNew->delta_F,stateNew->Fdot);

			//m_add(stateNew->F,stateOld->F,stateNew->delta_F);
			//sm_mlt(0.5,stateNew->delta_F,stateNew->delta_F);
			//m_inverse(stateNew->delta_F,stateNew->invF);

			//m_sub(stateNew->F, stateOld->F, stateNew->GRAD_U);
			//m_mlt(stateNew->GRAD_U,stateNew->invF,stateNew->h);


			// inverse deformation gradient			
			catchall(m_inverse(stateNew->F,stateNew->invF),
				printf("Inverse of F is singular for node %d F \n = %lf %lf %lf\n %lf %lf %lf\n %lf %lf %lf\n ",i,
					stateNew->F->me[0][0], stateNew->F->me[0][1], stateNew->F->me[0][2],
					stateNew->F->me[1][0], stateNew->F->me[1][1], stateNew->F->me[1][2],
					stateNew->F->me[2][0], stateNew->F->me[2][1], stateNew->F->me[2][2])
			);
	
			// // update Dn and Wn
			m_sub(stateNew->F, stateOld->F, stateNew->GRAD_U);

			m_mlt(stateNew->GRAD_U,stateOld->invF,stateNew->h);
			velocity_grad(stateNew->h, stateOld->L,stateOld->D, stateOld->W,dt,0);
			m_copy(stateOld->D,stateOld->Dbar);
	
			if ( dim == 2)
			{
				stateOld->div_v = stateOld->L->me[0][0] + stateOld->L->me[1][1];

			}else if ( dim == 3)
			{
				stateOld->div_v = stateOld->L->me[0][0] + stateOld->L->me[1][1] + stateOld->L->me[2][2];

			}

			if ( dim == 2)
			{
				stateOld->Dbar->me[0][0] += (-1.00/3.00)*stateOld->div_v;
				stateOld->Dbar->me[1][1] += (-1.00/3.00)*stateOld->div_v;

			}else {
				stateOld->Dbar->me[0][0] += (-1.00/3.00)*stateOld->div_v;
				stateOld->Dbar->me[1][1] += (-1.00/3.00)*stateOld->div_v;
				stateOld->Dbar->me[2][2] += (-1.00/3.00)*stateOld->div_v;

			}

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
			// m_inverse(stateNew->Fbar,stateNew->invFbar);
			// m_sub(stateNew->Fbar, stateOld->Fbar, stateNew->delta_F);

			// sm_mlt(1.000/dt,stateNew->delta_F,stateNew->Fbardot);
			// velocity_grad(stateNew->Lbar, stateNew->Dbar, stateNew->Wbar,stateNew->Fbardot,stateNew->invFbar);

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



			// find zi
			m_mlt(stateOld->V,stateOld->D,stateNew->temp);
			MAT * Vdot = m_get(dim,dim);
			double tempVec[dim];
			double omega[dim];
			double zi[dim];


			// zi = eigejk Vjm Dbamk
			for ( int i = 0  ; i < dim ; i++)
			{
				for ( int j = 0 ; j < dim ; j++)
				{
					for ( int k = 0 ; k < dim ; k++)
					{
						int e = ((i-j)*(j-k)*(k-i))/2;
						zi[i] += e*stateNew->temp->me[j][k];


					}
				}

			}


			m_ident(stateNew->temp);

			double traceV = stateOld->V->me[0][0] + stateOld->V->me[1][1] + stateOld->V->me[2][2];
			sm_mlt(traceV,stateNew->temp,stateNew->temp);
			m_sub(stateOld->V,stateNew->temp, stateNew->temp);
			m_inverse(stateNew->temp,stateNew->temp1);


			for ( int i = 0  ; i < dim ; i++)
			{
				for ( int j = 0 ; j < dim ; j++)
				{
					for ( int k = 0 ; k < dim ; k++)
					{
						int e = ((i-j)*(j-k)*(k-i))/2;
						omega[i] += e*stateOld->W->me[j][k];


					}
				}

			}

			for ( int i = 0  ; i < dim ; i++)
			{
				for ( int j = 0 ; j < dim ; j++)
				{
					tempVec[i] +=  stateNew->temp1->me[i][j]*zi[j]; 
				}

			}


			for ( int i = 0  ; i < dim ; i++)
			{
				omega[i] = omega[i] - 2*tempVec[i];		
			}

			for ( int i = 0  ; i < dim ; i++)
			{
				for ( int j = 0 ; j < dim ; j++)
				{
					for ( int k = 0 ; k < dim ; k++)
					{
						int e = ((i-j)*(j-k)*(k-i))/2;
						stateNew->Omega->me[i][j] += 0.5*e*omega[k];


					}
				}

			}

			m_mlt(stateOld->L,stateOld->V,Vdot);
			m_mlt(stateOld->V,stateNew->Omega,stateNew->temp);
			m_sub(Vdot,stateNew->temp,Vdot);

			sm_mlt(dt,Vdot,stateNew->V);
			m_add(stateNew->V,stateOld->V,stateNew->V);



			// Find rotation tensor using Hughes, Winget algorithm
			m_ident(stateNew->temp);
			sm_mlt(0.5*dt,stateNew->Omega,stateNew->temp1);
			m_sub(stateNew->temp,stateNew->temp1,stateNew->temp1);
			m_inverse(stateNew->temp1,stateNew->temp);

			m_ident(stateOld->temp);
			sm_mlt(0.5*dt,stateNew->Omega,stateNew->temp1);
			m_add(stateOld->temp,stateNew->temp1,stateNew->temp1);

			m_mlt(stateNew->temp1,stateOld->R,stateOld->temp);
			m_mlt(stateNew->temp,stateOld->temp,stateNew->R);
			// update rotation tensor
			printf("F = ");
			m_foutput(stdout,stateNew->F);
			printf("R = ");
			m_foutput(stdout,stateNew->R);
			printf("Vdot = ");
			m_foutput(stdout,Vdot);

			printf("V = ");
			m_foutput(stdout,stateNew->V);
			printf("W = ");
			m_foutput(stdout,stateOld->W);
			//printf("Vdot = ");
			//m_foutput(stdout,Vdot);
		







			// // Polar Decomposition
			// poldec(stateNew->Fbar, stateNew->R, stateNew->U, stateNew->V);

			//Get eigen values of polar decomposition
			// tracecatch(
			// symmeig(stateNew->V,stateNew->eigVecVBar,stateNew->eigValVBar);,
			// "Eigen values of V in internalForce");

			// tracecatch(
			// symmeig(stateNew->Dbar,stateNew->eigVecDBar,stateNew->eigValDBar);,
			// "Eigen values of D in internalForce");
			

			// // update critical network stretch 
			stateNew->critLambdaBar =lambdaCrit(stateOld->critLambdaBar,stateNew->eigValVBar, 
			stateOld->Dbar, critLambdaParams, stateNew->temperature, dt);

			stateNew->critLambdaBar = 2.4;



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
			stateOld->Jacobian = stateNew->Jacobian;
			stateOld->mSigma = stateNew->mSigma;			
			stateOld->critLambdaBar = stateNew->critLambdaBar;			
			m_copy(stateNew->V,stateOld->V);
			m_copy(stateNew->R,stateOld->R);






	return 0;
}