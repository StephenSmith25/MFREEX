#include "Deformation/update_Polar_Decomposition.h"


int update_Polar_Decomposition(state_variables * stateNew, state_variables * stateOld, double dt)
{



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
