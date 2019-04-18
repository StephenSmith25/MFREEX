#include "Material/material.h"




state_variables ** new_material_states(double * temperatures, int num_Points, int is_buckley,
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
		state[i]->delta_R = m_get(dim_s,dim_s);

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
		state[i]->sigma_R = m_get(dim_s,dim_s);



		// Material temperature
		if ( temperatures != NULL)
		{
			state[i]->temperature = temperatures[i]+273.15;
		}


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
			state[i]->dbar = m_get(dim_s,dim_s);

			// Conformation network stretch
			state[i]->Bbar = m_get(dim_s,dim_s);
			m_ident(state[i]->Bbar);
			state[i]->lambdaNBar = v_get(dim_s);
			state[i]->ep_n = m_get(dim_s,dim_s);


			state[i]->EP_bar = m_get(dim_s,dim_s);


			v_ones(state[i]->lambdaNBar);
			state[i]->Dn = m_get(dim_s,dim_s);
			state[i]->Ubar = m_get(dim_s,dim_s);
			m_ident(state[i]->Ubar);

			state[i]->delta_Ubar = m_get(dim_s,dim_s);

			state[i]->delta_ep_bar = m_get(dim_s,dim_s);

			state[i]->Vbar = m_get(dim_s,dim_s);
			m_ident(state[i]->Vbar);

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
		if ( is_plastic == 1)
		{
			state[i]->S_trial = m_get(dim_s,dim_s);
			state[i]->d_pl = m_get(dim_s,dim_s);
			state[i]->d_el = m_get(dim_s,dim_s);

			state[i]->gamma = 0;
			state[i]->Deps = 0;
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


state_variables * new_material_state( double temperature, MATERIAL_TYPE mat_type, int dim, int is_AXI)

{
		int dim_s = dim;
		if ( is_AXI == 1)
		{
		dim_s = 3;
		}

		state_variables * state =  malloc(sizeof(state_variables));

		/* ------------------------------------------*/
		/* -----------------Deformation--------------*/
		/* ------------------------------------------*/


		// Deformation gradient
		state->F = m_get(dim_s,dim_s);
		m_ident(state->F);
		state->invF = m_get(dim_s,dim_s);
		m_ident(state->invF);
		state->Jacobian = 1.00;
		state->Fdot = m_get(dim_s,dim_s);

		// Deformation tensors
		state->B = m_get(dim_s,dim_s);
		state->C = m_get(dim_s,dim_s);

		// Polar Decomposition
		state->R = m_get(dim_s,dim_s);
		state->delta_R = m_get(dim_s,dim_s);

		m_ident(state->R);
		state->U = m_get(dim_s,dim_s);
		state->V = m_get(dim_s,dim_s);
		m_ident(state->U);
		m_ident(state->V);
		state->Vdot = m_get(dim_s,dim_s);


		// Rate measures	
		state->L = m_get(dim_s,dim_s);
		state->D = m_get(dim_s,dim_s);
		state->W = m_get(dim_s,dim_s);
		state->Omega = m_get(dim_s,dim_s);
		state->div_v = 0;



		// Rotated measures
		state->d = m_get(dim_s,dim_s);
		state->omega = m_get(dim_s,dim_s);
		state->T = m_get(dim_s,dim_s);



		// Stress
		state->sigma = m_get(dim_s,dim_s);
		state->sigma_R = m_get(dim_s,dim_s);



		// Material temperature
		state->temperature = temperature+273.15;


		/* ------------------------------------------*/
		/* ------------------Buckley-----------------*/
		/* ------------------------------------------*/

		if (mat_type == BUCKLEY)
		{

			// Eigen values

			state->Fbar = m_get(dim_s,dim_s);
			m_ident(state->Fbar);
			state->Dbar = m_get(dim_s,dim_s);
			state->Lbar = m_get(dim_s,dim_s);
			state->dbar = m_get(dim_s,dim_s);

			// Conformation network stretch
			state->Bbar = m_get(dim_s,dim_s);
			m_ident(state->Bbar);
			state->lambdaNBar = v_get(dim_s);
			state->ep_n = m_get(dim_s,dim_s);


			state->EP_bar = m_get(dim_s,dim_s);


			v_ones(state->lambdaNBar);
			state->Dn = m_get(dim_s,dim_s);
			state->Ubar = m_get(dim_s,dim_s);
			m_ident(state->Ubar);

			state->delta_Ubar = m_get(dim_s,dim_s);

			state->delta_ep_bar = m_get(dim_s,dim_s);

			state->Vbar = m_get(dim_s,dim_s);
			m_ident(state->Vbar);

			// Eigen values

			state->eigVecDBar = m_get(dim_s,dim_s);
			state->eigValDBar = v_get(dim_s);
			state->eigVecVBar = m_get(dim_s,dim_s);
			state->eigValVBar = v_get(dim_s);
			state->lambdaDot = v_get(dim_s);


			// Buckley output stresses
			state->Sb = m_get(dim_s,dim_s);
			state->Sc = m_get(dim_s,dim_s);
			state->mSigma = 0;


			// rotated stresses
			state->Sb_R = m_get(dim_s,dim_s);
			state->Sc_R = m_get(dim_s,dim_s);


			// Material constants
			state->gamma = 0.653e6;
			state->critLambdaBar = 100.00;
			state->tau = 0;
			state->lambdaNMax = 1;



		}

		/* ------------------------------------------*/
		/* ---------------Plasticity-----------------*/
		/* ------------------------------------------*/
		if ( mat_type == PLASTIC)
		{
			state->S_trial = m_get(dim_s,dim_s);
			state->d_pl = m_get(dim_s,dim_s);
			state->d_el = m_get(dim_s,dim_s);

			state->gamma = 0;
			state->Deps = 0;
			// do nothing for now
		}

		/* ------------------------------------------*/
		/* -------------Workspace arrays-------------*/
		/* ------------------------------------------*/

		state->v_temp1 = v_get(dim_s);
		state->v_temp2 = v_get(dim_s);
		state->v_temp3 = v_get(dim_s);
		state->v_temp4 = v_get(dim_s);


		state->m_temp1 = m_get(dim_s,dim_s);
		state->m_temp2 = m_get(dim_s,dim_s);
		state->m_temp3 = m_get(dim_s,dim_s);
		state->m_temp4 = m_get(dim_s,dim_s);



	


		return state;



}




MATERIAL * create_new_material(MATERIAL_TYPE material_type, void * MATERIAL_LAW, VEC * params)
{

	MATERIAL * material = malloc(1*sizeof(MATERIAL));

	switch (material_type)
	{
		case (HYPERELASTIC):
		{
			HYPERELASTIC_LAW hyperelastic_law = *((HYPERELASTIC_LAW*) MATERIAL_LAW);

			switch(hyperelastic_law)
			{
				case (MOONEY_RIVLIN):
				{
					MOONEY_RIVLIN_MATERIAL * material_rivlin = malloc(1*sizeof(MOONEY_RIVLIN_MATERIAL));
					//material_rivlin->params = params;
					material->MATERIAL_LAW = material;
					printf("MOONEY_RIVLIN LAW");
					break;
				}
				case (NEO_HOOKEAN):
				{
					break;
				} 
				default:
				{
					fprintf(stderr,"Material not yet supported \n");
				}
			}

			break;
		}
		case (PLASTIC):
		{


			break;
		}
		case ( BUCKLEY):
		{


			break;
		}
		default:
		{
			fprintf(stderr,"material model not supported\n");
			return NULL;
		}
	}


	return material;
}
