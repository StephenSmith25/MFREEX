#include "Material/material.h"




state_variables ** new_material_state(double * temperatures, int num_Points, int is_buckley,
	int is_plastic, int dim, int is_AXI)
{


	state_variables ** state = malloc(num_Points * sizeof(state_variables*));

	int dim_s = dim;
	if ( is_AXI == 1)
	{
		dim_s = 3;
	}

	

	for ( int i = 0 ; i < num_Points; i++){

		state[i] =  malloc(sizeof(state_variables));

		/* ------------------------------------------*/
		/* -----------------Deformation--------------*/
		/* ------------------------------------------*/


		// Deformation gradient
		state[i]->F = m_get(dim_s,dim_s);
		m_ident(state[i]->F);
		state[i]->invF = m_get(dim_s,dim_s);
		m_ident(state[i]->invF);
		state[i]->Jacobian = 1.00;

		// Deformation tensors
		state[i]->B = m_get(dim_s,dim_s);
		state[i]->C = m_get(dim_s,dim_s);

		// Polar Decomposition
		state[i]->R = m_get(dim_s,dim_s);
		m_ident(state[i]->R);
		state[i]->U = m_get(dim_s,dim_s);
		state[i]->V = m_get(dim_s,dim_s);
		m_ident(state[i]->U);
		m_ident(state[i]->V);
		state[i]->Vdot = m_get(dim_s,dim_s);


		// Rate measures	
		state[i]->L = m_get(dim_s,dim_s);
		state[i]->D = m_get(dim_s,dim_s);
		state[i]->W = m_get(dim_s,dim_s);
		state[i]->Omega = m_get(dim_s,dim_s);
		state[i]->div_v = 0;



		// Rotated measures
		state[i]->d = m_get(dim_s,dim_s);
		state[i]->omega = m_get(dim_s,dim_s);
		state[i]->T = m_get(dim_s,dim_s);



		// Stress
		state[i]->sigma = m_get(dim_s,dim_s);



		// Material temperature
		state[i]->temperature = temperatures[i]+273.15;



		/* ------------------------------------------*/
		/* ------------------Buckley-----------------*/
		/* ------------------------------------------*/

		if (is_buckley == 1)
		{

			// Eigen values

			state[i]->Fbar = m_get(dim_s,dim_s);
			m_ident(state[i]->Fbar);
			state[i]->Dbar = m_get(dim_s,dim_s);
			state[i]->Lbar = m_get(dim_s,dim_s);


			// Conformation network stretch
			state[i]->Bbar = m_get(dim_s,dim_s);
			m_ident(state[i]->Bbar);
			state[i]->lambdaNBar = v_get(dim_s);
			state[i]->Dn = m_get(dim_s,dim_s);


			// Eigen values

			state[i]->eigVecDBar = m_get(dim_s,dim_s);
			state[i]->eigValDBar = v_get(dim_s);
			state[i]->eigVecVBar = m_get(dim_s,dim_s);
			state[i]->eigValVBar = v_get(dim_s);
			state[i]->lambdaDot = v_get(dim_s);


			// Buckley output stresses
			state[i]->Sb = m_get(dim_s,dim_s);
			state[i]->Sc = m_get(dim_s,dim_s);
			state[i]->mSigma = 0;


			// rotated stresses
			state[i]->Sb_R = m_get(dim_s,dim_s);
			state[i]->Sc_R = m_get(dim_s,dim_s);


			// Material constants
			state[i]->gamma = 0.653e6;
			state[i]->critLambdaBar = 100.00;
			state[i]->tau = 0;
			state[i]->lambdaNMax = 1;



		}

		/* ------------------------------------------*/
		/* ---------------Plasticity-----------------*/
		/* ------------------------------------------*/
		if ( is_plastic == 0)
		{

			// do nothing for now
		}

		/* ------------------------------------------*/
		/* -------------Workspace arrays-------------*/
		/* ------------------------------------------*/

		state[i]->v_temp1 = v_get(dim_s);
		state[i]->v_temp2 = v_get(dim_s);
		state[i]->v_temp3 = v_get(dim_s);
		state[i]->v_temp4 = v_get(dim_s);


		state[i]->m_temp1 = m_get(dim_s,dim_s);
		state[i]->m_temp2 = m_get(dim_s,dim_s);
		state[i]->m_temp3 = m_get(dim_s,dim_s);
		state[i]->m_temp4 = m_get(dim_s,dim_s);



	}




	return state;
}
