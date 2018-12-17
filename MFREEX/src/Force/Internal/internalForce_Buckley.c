#include "Force/Internal/internalForce_Buckley.h"


static int call_count;

double internalForce_hyperelastic(VEC * Fint, SCNI_OBJ * scni_obj, VEC * disp, VEC * velocity,VEC * matParams, 
	state_Buckley ** state_n ,int is_axi, int dim, double deltat){


	__zero__(Fint->ve,Fint->max_dim);



	int dim_v = 0;
	int dim_s = 0;
	// check if problem is axisymmetric
	if ( is_axi == 1){
		dim_v = 5;
		dim_s = 3;
	}else{
		dim_v = dim*dim;
		dim_s = dim;
	}

	// loop over all integration points

	SCNI ** scni = scni_obj->scni;
	int num_int_points = scni_obj->num_points;

	
	// time step calculation
	double delta_t_min = 1000;


	// Find effective lambda and mu
	double lambda = matParams->ve[0];
	double mu = matParams->ve[1];

	// set number of threads
	omp_set_num_threads(8);


#pragma omp parallel 
	{

		MAT * F_n_1 = m_get(dim_s,dim_s);
		MAT * F_n_1_bar = m_get(dim_s,dim_s);
		MAT * invF_n_1 = m_get(dim_s,dim_s);
		MAT * invF_n_1_bar = m_get(dim_s,dim_s);
		MAT * Fdot = m_get(dim_s,dim_s);
		MAT * Fdot_bar = m_get(dim_s,dim_s);
		VEC * stressVoigt = v_get(dim_v);
		VEC * fIntTemp = v_get(Fint->max_dim);
		MAT * R = m_get(dim_s,dim_s);
		MAT * U = m_get(dim_s,dim_s);
		MAT * V = m_get(dim_s,dim_s);	
		MAT * temp = m_get(dim_s,dim_s);
		MAT * delta_F = m_get(dim_s,dim_s);


		// Stress variables
		MAT * cauchy_dev = m_get(dim_s,dim_s);
		MAT * cauchy_hyd = m_get(dim_s,dim_s);

		// strain varaibles
		MAT * delta_ep_true = m_get(dim_s,dim_s);
		MAT * delta_ep_true_vol = m_get(dim_s,dim_s);
		MAT * delta_ep_true_iso = m_get(dim_s,dim_s);


		MAT * delta_cauchy_dev = m_get(dim_s,dim_s);
		MAT * delta_cauchy_hyd = m_get(dim_s,dim_s);
		MAT * Sb_n_1 = m_get(dim_s,dim_s);
		MAT * Sc_n_1 = m_get(dim_s,dim_s);

		double div_v;
		double delta_t_min_i = 1000;
		MAT * B;
		IVEC * neighbours;
		MAT * F_r;
		double Jacobian_n;
		double Jacobian_n_1;

#pragma omp for nowait schedule(dynamic,4) 
		for(int i = 0 ; i < num_int_points ; i++){



			/* ------------------------------------------*/
			/* -----------------Deformation--------------*/
			/* ------------------------------------------*/

			//------------------------//
			//        Def Grad        //
			//------------------------//


			/*  Find deformation gradient */
			B = scni[i]->B;
			neighbours = scni[i]->sfIndex;
			int num_neighbours = neighbours->max_dim;
			
			/*  Find deformation gradient */
			get_defgrad(F_n_1, B, neighbours,F_r, disp);

			/* Find Fdot and Fbar dot */
			m_sub(F_n_1, state_n[i]->F, delta_F);
			__smlt__(delta_F->base, 1.00/deltat, Fdot->base, dim_s*dim_s);

			// inverse deformation gradient
			m_inverse(F_n_1,invF_n_1);

			// polar decomposition to find delta R
			poldec(delta_F, R, U, V);

			//------------------------//
			//      Velocity grad     //
			//------------------------//

			// Find velocity gradient 
			velocity_grad(state_n[i]->L, state_n[i]->D, state_n[i]->W, Fdot,invF_n_1);

			// div_v

			if ( dim_s == 2)
			{
				div_v = state_n[i]->L->me[0][0] + state_n[i]->L->me[1][1];

			}else if ( dim_s == 3)
			{
				div_v = state_n[i]->L->me[0][0] + state_n[i]->L->me[1][1] + state_n[i]->L->me[2][2];

			}

			state_n[i]->div_v = div_v;
			//------------------------//
			//   Update true strain   //
			//------------------------//


			sm_mlt(deltat, state_n[i]->D, delta_ep_true);
			// update total strain
			m_mlt(R, state_n[i]->ep_true, temp);
			mmtr_mlt(temp, R, state_n[i]->ep_true);
			m_add(state_n[i]->ep_true, delta_ep_true, state_n[i]->ep_true);
			double delta_ep_vol = 0;

			if ( dim == 2)
			{
				delta_ep_vol = (1.00/3.00)*(delta_ep_true->me[0][0] +delta_ep_true->me[1][1]);
				delta_ep_true_vol->me[0][0] = delta_ep_vol;
				delta_ep_true_vol->me[1][1] = delta_ep_vol;

			}else if ( dim == 3)
			{
				delta_ep_vol = (1.00/3.00)*(delta_ep_true->me[0][0] +delta_ep_true->me[1][1] +delta_ep_true->me[2][2] );
				delta_ep_true_vol->me[0][0] = delta_ep_vol;
				delta_ep_true_vol->me[1][1] = delta_ep_vol;
				delta_ep_true_vol->me[2][2] = delta_ep_vol;
			}

			m_sub(delta_ep_true,delta_ep_true_vol,delta_ep_true_iso);



			// Update Jacobian
			Jacobian_n = state_n[i]->Jacobian;
			Jacobian_n_1 = Jacobian_n + Jacobian_n*div_v*deltat;

			/* ------------------------------------------*/
			/* ---------Isochoric   Deformation---------*/
			/* ------------------------------------------*/

			/* Distortional deformation */
			__smlt__(F_n_1->base, pow(Jacobian_n_1,-1.00/3.00), F_n_1_bar->base, dim_s*dim_s);
			// inv Fbar_n_1
			__smlt__(invF_n_1->base,pow(Jacobian_n_1,1.00/3.00), invF_n_1_bar->base, dim_s*dim_s);

			// F bar dot
			m_sub(F_n_1_bar, state_n[i]->Fbar, Fdot_bar);
			__smlt__(Fdot->base, 1.00/deltat, Fdot_bar->base, dim_s*dim_s);

			velocity_grad(state_n[i]->Lbar, state_n[i]->Dbar, state_n[i]->Wbar, Fdot_bar, invF_n_1_bar);


			/* ------------------------------------------*/
			/* ------------- --Find Stress--------------*/
			/* ------------------------------------------*/
			/*  Obtain stresses using explicit integration of stress rate */
			buckleyBond(Sb_n_1,state_n[i],matParams,deltat);

			// Conformational stress
			buckleyConf(Sc_n_1,state_n[i],matParams,deltat,Jacobian_n_1);



			/* ------------------------------------------*/
			/* --------------New time increment----------*/
			/* ------------------------------------------*/
			m_add(Sb_n_1, Sc_n_1, cauchy_dev);


			// Find hydrostatic and deviatoric stress

			m_sub(cauchy_dev,state_n[i]->dev_Stress,delta_cauchy_dev);

			for ( int k = 0 ; k < dim_s ; k++)
			{
				cauchy_hyd->me[i][i] = matParams->ve[5]*log(Jacobian_n_1);

			}
			m_sub(cauchy_hyd,state_n[i]->hyd_Stress,delta_cauchy_hyd);


			/* ------------------------------------------*/
			/* ---------------Bulk Damping---------------*/
			/* ------------------------------------------*/



			/* ------------------------------------------*/
			/* --------------Internal Force--------------*/
			/* ------------------------------------------*/


			/* Integration parameter */
			double intFactor = scni[i]->area;
			if( is_axi == 1){

				intFactor = intFactor * 2*PI*scni[i]->r;
			}

			// Find bond and confromation stresss

			// Find increments 
			// Find Delta_



			//Find 1D frequency bounds
			// time step calculation
			double B11 = 0;
			double B22 = 0;
			double MaxB = -1;
			double rho = 1000e-9;


			for ( int k = 0 ; k < (num_neighbours)*2 ; k++)
			{	
				B11 += (num_neighbours/rho)*(B->me[0][k]*B->me[0][k]);
				B22 += (num_neighbours/rho)*(B->me[1][k]*B->me[1][k]);

			}

			MaxB = max(B11,B22);



			double b1 = 0.08;
			double b2 = 1.44;
			double Le = 1.6;
			double c = sqrt(((lambda+2*mu)/rho));

			double P_b1 = b1*div_v*rho*Le*c;
			double eta = b1;
			double P_b2 = 0;
			if ( div_v < 0 ){
				P_b2 = Le*rho*(b2*b2*Le*div_v*div_v);
				eta -= b2*b2*Le*(1/c)*div_v;
			}


			double delta_t = (2.00/sqrt(((lambda + 2*mu)*MaxB)))*(sqrt(1+eta*eta)-eta);
			if ( delta_t <delta_t_min_i )
			{
				delta_t_min_i = delta_t;
			}



			__zero__(scni[i]->fInt->ve,scni[i]->fInt->max_dim);


			// // push forward stress to new reference configuration
			// if ( dim_v == 4)
			// {


			// 	double S11_hyd = (Jacobian)*(P_b1+P_b2)*invF11;
			// 	double S22_hyd = (Jacobian)*(P_b1+P_b2)*invF22;


			// 	S11 = stressVoigt->ve[0]+S11_hyd; S12 = stressVoigt->ve[2];
			// 	S21 = stressVoigt->ve[3];S22 = stressVoigt->ve[1]+S22_hyd;

			// 	F11 = scni[i]->F_r->me[0][0]; F12 = scni[i]->F_r->me[0][1];
			// 	F21 = scni[i]->F_r->me[1][0]; F22 = scni[i]->F_r->me[1][1];

			// 	stressVoigt->ve[0] = S11*F11 + S12*F12;
			// 	stressVoigt->ve[1] = S21*F21 + S22*F22; 
			// 	stressVoigt->ve[2] = S11*F21 + S12*F22;
			// 	stressVoigt->ve[3] = S21*F11 + S22*F12;

			// }else if ( dim_v == 5)
			// {
			// 	S11 = stressVoigt->ve[0]; S12 = stressVoigt->ve[2];
			// 	S21 = stressVoigt->ve[3]; S22 = stressVoigt->ve[1];
			// 	S33 = stressVoigt->ve[5];

			// 	F11 = scni[i]->F_r->me[0][0]; F12 = scni[i]->F_r->me[0][1];
			// 	F21 = scni[i]->F_r->me[1][0]; F22 = scni[i]->F_r->me[1][1];
			// 	F33 = scni[i]->F_r->me[2][2];

			// 	stressVoigt->ve[0] = S11*F11 + S12*F12;
			// 	stressVoigt->ve[1] = S21*F21 + S22*F22; 
			// 	stressVoigt->ve[2] = S11*F21 + S12*F22;
			// 	stressVoigt->ve[3] = S21*F11 + S22*F12;
			// 	stressVoigt->ve[4] = S33*F33;

			// }else{
			// 	fprintf(stderr,"dimension not yet supported");
			// }


			vm_mlt(scni[i]->B,stressVoigt,scni[i]->fInt);

			for ( int k = 0 ; k < num_neighbours; k++){
				fIntTemp->ve[2*scni[i]->sfIndex->ive[k]] += intFactor * scni[i]->fInt->ve[2*k];
				fIntTemp->ve[2*scni[i]->sfIndex->ive[k]+1] += intFactor * scni[i]->fInt->ve[2*k+1];
			}
		}

	/*  Make this atomic or mutex so it is only done by one thread */
	#pragma omp critical
		{
			__add__(fIntTemp->ve, Fint->ve, Fint->ve, Fint->max_dim);
			if ( delta_t_min_i < delta_t_min)
			{
				delta_t_min = delta_t_min_i;
			}
		}

	// /*  Free allocated memory */
	// 	v_free(stressVoigt);
	// 	v_free(fIntTemp);
	// 	M_FREE(F);
	// 	M_FREE(Fdot);

	}
	++call_count;

	return delta_t_min;
}
