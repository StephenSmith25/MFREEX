#include "Deformation/update_Polar_Decomposition.h"
#include "dsyevq3.h"
#include "dsyevv3.h"
#include "dsyevh3.h"
#include <math.h>


int update_Polar_Decomposition(state_variables * stateNew, state_variables * stateOld, double dt)
{


			// algorithm from simo
			// mmtr_mlt(stateNew->F, stateNew->F, stateNew->B);
			// MAT * Q = stateNew->m_temp1;
			// double w[3];
			// dsyevh3(stateNew->B->me,Q->me,w);

			// MAT * V_n_1 = stateNew->V;
			// MAT * V_n = stateOld->V;
			// MAT * B = stateNew->B;
			// MAT * invV_n_1 = stateNew->m_temp4;

			// double lambda_1 = sqrt(w[0]);
			// double lambda_2 = sqrt(w[1]);
			// double lambda_3 = sqrt(w[2]);

			// double i1 = lambda_1 + lambda_2 + lambda_3;
			// double i2 = lambda_1*lambda_2 + lambda_1*lambda_3 + lambda_2*lambda_3;
			// double i3 = lambda_1*lambda_2*lambda_3;

			// double D = i1*i2 - i3;
			// double invD = 1.000/D;


			// // FInd Bsqr 

			// MAT * B_sqr = stateNew->m_temp2;

			// m_zero(B_sqr);

			// B_sqr->me[0][0] = lambda_1*lambda_1;
			// B_sqr->me[1][1] = lambda_2*lambda_2;
			// B_sqr->me[2][2] = lambda_3*lambda_3;





			// // Find V_n_1 = 1/D [ -B^2 + (i1^2 - i2)B + i1*i3 *1]
			// double fac1 = (i1*i1 - i2);
			// double fac2 = i1*i3;
			// double fac3 = 1.00/i3;

			// for ( int i = 0 ; i < 3 ; i++)
			// {
			// 	for ( int j = 0 ; i < 3 ; j++)
			// 	{
			// 		if ( i == j)
			// 		{
			// 			V_n_1->me[i][j] = invD*(-B_sqr->me[i][j] + fac1*B->me[i][j] + fac2 );
			// 			invV_n_1->me[i][j] = fac3*(B->me[i][j] - i1*V_n_1->me[i][j] + i2);

			// 		}else{
			// 			V_n_1->me[i][j] = invD*(-B_sqr->me[i][j] + fac1*B->me[i][j] );
			// 			invV_n_1->me[i][j] = fac3*(B->me[i][j] - i1*V_n_1->me[i][j] );

			// 		}
			// 	}
			// }
			// // Find V-1
			// stateNew->
			// // FInd R




			// FIND THE ROTATION TENSOR
			const int dim = 3;
			//int dim = stateOld->V->m;
			m_mlt(stateOld->V,stateNew->D,stateNew->m_temp1);
			m_mlt(stateNew->D,stateOld->V,stateNew->m_temp2);
			m_sub(stateNew->m_temp2,stateNew->m_temp1,stateNew->m_temp1);

			MAT * Vdot = stateNew->Vdot;
			VEC * h = stateNew->v_temp1;
			VEC * w = stateNew->v_temp2;
			VEC * omega = stateNew->v_temp3;
			VEC * z = stateNew->v_temp4;

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
						z->ve[j] += (double)(0.5*e)*stateNew->m_temp1->me[i][k];
						w->ve[j] += (double)(0.5*e)*stateNew->W->me[i][k];

					}
				}

			}


			m_ident(stateNew->m_temp1);

			double traceV = stateOld->V->me[0][0] + stateOld->V->me[1][1] + stateOld->V->me[2][2];
			sm_mlt(traceV,stateNew->m_temp1,stateNew->m_temp1);
			m_sub(stateNew->m_temp1,stateOld->V, stateNew->m_temp1);

			// m_inverse_small(stateNew->m_temp1,stateNew->m_temp4);
			// vm_mlt(stateNew->m_temp4, z, h);

			CHfactor(stateNew->m_temp1);
			CHsolve(stateNew->m_temp1, z, h);
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
			m_mlt(stateOld->V,stateNew->Omega,stateNew->m_temp1);
			m_sub(Vdot,stateNew->m_temp1,Vdot);


			sm_mlt(dt,Vdot,stateNew->V);
			m_add(stateNew->V,stateOld->V,stateNew->V);

			// //Find rotation tensor using Hughes, Winget algorithm
			m_ident(stateNew->m_temp1);
			sm_mlt(0.5*dt,stateNew->Omega,stateNew->m_temp2);
			m_sub(stateNew->m_temp1,stateNew->m_temp2,stateNew->m_temp2);
			m_inverse_small(stateNew->m_temp2,stateNew->m_temp1);

			m_ident(stateOld->m_temp1);
			sm_mlt(0.5*dt,stateNew->Omega,stateNew->m_temp2);
			m_add(stateOld->m_temp1,stateNew->m_temp2,stateNew->m_temp2);

			m_mlt(stateNew->m_temp2,stateOld->R,stateOld->m_temp1);
			m_mlt(stateNew->m_temp1,stateOld->m_temp1,stateNew->R);








	return 0.0;
}
int polar_decomposition(state_variables * stateNew, state_variables * stateOld, double dt)
{






	return 0;
}
