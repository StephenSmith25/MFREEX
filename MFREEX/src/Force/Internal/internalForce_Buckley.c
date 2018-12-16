#include "Force/Internal/internalForce_Buckley.h"




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
		MAT * invF_n_h = m_get(dim_s,dim_s);
		MAT * invF_n_1 = m_get(dim_s,dim_s);
		MAT * Fdot_n_h = m_get(dim_s,dim_s);
		MAT * F_n_h = m_get(dim_s,dim_s);

		VEC * stressVoigt = v_get(dim_v);
		VEC * fIntTemp = v_get(Fint->max_dim);
		
		double div_v;
		double delta_t_min_i = 1000;
		MAT * B ;
		MAT * F_n;
		MAT * L_n_h;
		MAT * D_n_h;
		MAT * W_n_h;
		IVEC * neighbours;
		MAT * F_r;
		double Jacobian_n;
		double Jacobian_n_1;

#pragma omp for nowait schedule(dynamic,4) 
		for(int i = 0 ; i < num_int_points ; i++){

			/*  Find deformation gradient */
			B = scni[i]->B;
			neighbours = scni[i]->sfIndex;
			int num_neighbours = neighbours->max_dim;
			
			/*  Find deformation gradient */
			get_defgrad(F_n_1, B, neighbours,F_r, disp);
			get_dot_defgrad(Fdot_n_h,B, neighbours,F_r,velocity);

			// inverse deformation gradient
			m_inverse(F_n_1,invF_n_1);
			F_n = state_n[i]->F;

			// Find half time step deformation gradient 
			ms_mltadd(F_n, Fdot_n_h, 0.5*deltat, F_n_h);
			m_inverse(F_n_h,invF_n_h);
			
			L_n_h = state_n[i]->L;
			D_n_h = state_n[i]->D;
			W_n_h = state_n[i]->W;

			// Find velocity gradient 
			velocity_grad(L_n_h, D_n_h, W_n_h, Fdot_n_h,invF_n_h);

			// Find Jacobian at t_n_1;

			if ( dim_s == 2)
			{
				div_v = L_n_h->me[0][0] + L_n_h->me[1][1];

			}else if ( dim_s == 3)
			{
				div_v = L_n_h->me[0][0] + L_n_h->me[1][1] + L_n_h->me[2][2];

			}
			// Update Jacobian
			Jacobian_n = state_n[i]->Jacobian;
			Jacobian_n_1 = Jacobian_n + Jacobian_n*div_v*deltat;

	
			/* Integration parameter */
			double intFactor = scni[i]->area;
			if( is_axi == 1){

				intFactor = intFactor * 2*PI*scni[i]->r;
			}

			// Find bond and confromation stresss










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

	/*  Free allocated memory */
		v_free(stressVoigt);
		v_free(fIntTemp);
		M_FREE(F);
		M_FREE(Fdot);

	}
	++call_count;

	return delta_t_min;
}
