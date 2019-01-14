#include "Force/Internal/internalForce_hyperelastic.h"

static int call_count ;

double internalForce_hyperelastic(VEC * Fint, SCNI_OBJ * scni_obj, VEC * disp, VEC * velocity, VEC * matParams, char * Material, int is_axi, int dim){




	__zero__(Fint->ve,Fint->max_dim);

	double lambda, mu;
	// create pointer to correct material
	int (*mat_func_ptr)(VEC*,MAT*,VEC*) = NULL;

	if ( strcmp (Material, "SVK") == 0 )
	{
		mat_func_ptr = &SVK;
		lambda = matParams->ve[0];
		mu = matParams->ve[1];

	}else if( strcmp(Material, "YEOH") == 0 )
	{
		mat_func_ptr = &yeoh;
		mu = 2*(matParams->ve[0]);
		lambda = 2*matParams->ve[3] - (2/3)*mu;

	}

	int dim_v = 0;
	// check if problem is axisymmetric
	if ( is_axi == 1){
		dim_v = 5;
		dim = 3;
	}else{
		dim_v = dim*dim;
	}

	// loop over all integration points

	SCNI ** scni = scni_obj->scni;
	int num_int_points = scni_obj->num_points;

	
	// time step calculation
	// 
	double delta_t_min = 1000;



	omp_set_num_threads(3);



#pragma omp parallel 
{

	MAT * F = m_get(dim,dim);
	MAT * Fdot = m_get(dim,dim);
	MAT * L = m_get(dim,dim);
	VEC * stressVoigt = v_get(dim_v);
	MAT * invF = m_get(dim,dim);
	VEC * fIntTemp = v_get(Fint->max_dim);
	// stress tensor
	double S11,S12,S13,S21,S22,S23,S31,S32,S33;
	double F11,F12,F13,F21,F22,F23,F31,F32,F33;
	double invF11,invF12,invF13,invF21,invF22,invF23,invF31,invF32,invF33;
	double Co11,Co12,Co13,Co21,Co22,Co23,Co31,Co32,Co33;
	double L11,L12,L13,L21,L22,L23,L31,L32,L33;
	double div_v;
	double delta_t_min_i = 1000;
	MAT * B ;
	IVEC * neighbours;
	MAT * F_r;


#pragma omp for nowait schedule(dynamic,4) 
	for(int i = 0 ; i < num_int_points ; i++){

		/*  Find deformation gradient */
		B = scni[i]->B;
		neighbours = scni[i]->sfIndex;
		F_r = scni[i]->F_r;
		/*  Find deformation gradient */
		get_defgrad(F, B, neighbours,F_r, disp);
		get_dot_defgrad(Fdot,B, neighbours,F_r,velocity);

		int num_neighbours = neighbours->max_dim;
	




		double Jacobian ;


		if ( dim == 2)
		{	

			F11 = F->me[0][0]; F12 = F->me[0][1];
			F21 = F->me[1][0]; F22 = F->me[1][1];


			Jacobian = F11*F22 - F12*F21;

			Co11 = F22;
			Co12 = -F21;
			Co21 = -F12;
			Co22 = F11;
			invF11 = Co11/Jacobian; invF12 = Co21/Jacobian;
			invF21 = Co12/Jacobian; invF22 = Co22/Jacobian;


			// Find velocity gradient

			L11 = Fdot->me[0][0]*invF11 + Fdot->me[0][1]*invF21; L12 = Fdot->me[0][0]*invF12 + Fdot->me[0][1]*invF22;
			L21 = Fdot->me[1][0]*invF11 + Fdot->me[1][1]*invF21; L22 = Fdot->me[1][0]*invF12 + Fdot->me[1][1]*invF22;



			div_v = L11 + L22;
		}else if ( dim == 3)
		{

			F11 = F->me[0][0]; F12 = F->me[0][1]; F12 = F->me[0][1];
			F21 = F->me[1][0]; F22 = F->me[1][1]; F12 = F->me[0][1];
			F21 = F->me[1][0]; F22 = F->me[1][1]; F12 = F->me[0][1];
			m_inverse(F,invF);

			m_mlt(Fdot,invF,L);

			Jacobian = determinant(F);

			// Find velocity gradient
			div_v = L->me[0][0] + L->me[1][1] + L->me[2][2];	
			div_v = L->me[0][0] + L->me[1][1];
		}
		/* Integration parameter */
		double intFactor = scni[i]->area;
		if( is_axi == 1){

			intFactor = intFactor * 2*PI*scni[i]->r;
		}

		//Find 1D frequency bounds

		double B11 = 0;
		double B22 = 0;
		double B33 = 0;

		double MaxB = -1;
		double rho = 1380e-9;
		double num_neighbours_J = 0;
		double m_j;
		double volume_I = intFactor;
		for ( int k = 0 ; k < (num_neighbours) ; k++)
		{	
			int index_J = neighbours->ive[k];
			num_neighbours_J = scni[index_J]->sfIndex->max_dim;
			volume_I = scni[i]->area;
			m_j = rho * scni[index_J]->area;

			// if ( is_axi == 1){
			// 	m_j = rho * scni[index_J]->area*2*PI*scni[index_J]->r;
			// }

			B11 += (volume_I*num_neighbours_J/m_j)*(B->me[0][2*k]*B->me[0][2*k]);
			B22 += (volume_I*num_neighbours_J/m_j)*(B->me[1][2*k+1]*B->me[1][2*k+1]);

			if ( is_axi == 1)
			{
				B33 += (volume_I*num_neighbours_J/m_j)*(B->me[4][2*k]*B->me[4][2*k]);

			}

		}

		MaxB = max(B11,B22);

		// if ( is_axi == 1)
		// {
		// 	MaxB = max(MaxB,B33);

		// }


		double b1 = 0.06;
		double b2 = 1.44;

		double Le = sqrt(scni[i]->area);
		//Le = 3.5;
		double c = sqrt(((lambda+2*mu)/rho));

		double P_b1 = 0;
		P_b1 = b1*div_v*rho*Le*c;

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



		mat_func_ptr(stressVoigt,F,matParams);
		__zero__(scni[i]->fInt->ve,scni[i]->fInt->max_dim);


		// push forward stress to new reference configuration
		if ( dim_v == 4)
		{


			double S11_hyd = (Jacobian)*(P_b1+P_b2)*invF11;
			double S22_hyd = (Jacobian)*(P_b1+P_b2)*invF22;


			S11 = stressVoigt->ve[0]+S11_hyd; S12 = stressVoigt->ve[2];
			S21 = stressVoigt->ve[3]; S22 = stressVoigt->ve[1]+S22_hyd;

			F11 = scni[i]->F_r->me[0][0]; F12 = scni[i]->F_r->me[0][1];
			F21 = scni[i]->F_r->me[1][0]; F22 = scni[i]->F_r->me[1][1];

			stressVoigt->ve[0] = S11*F11 + S12*F12;
			stressVoigt->ve[1] = S21*F21 + S22*F22; 
			stressVoigt->ve[2] = S11*F21 + S12*F22;
			stressVoigt->ve[3] = S21*F11 + S22*F12;

		}else if ( dim_v == 5)
		{
			double S11_hyd = (Jacobian)*(P_b1+P_b2)*invF->me[0][0];
			double S22_hyd = (Jacobian)*(P_b1+P_b2)*invF->me[1][1];
			double S33_hyd = (Jacobian)*(P_b1+P_b2)*invF->me[2][2];
			S33_hyd = 0;

			if (( i == 0) && (call_count % 1000 == 0))
			{
				printf("S11_HYD = %lf \n", S11_hyd);
				printf("S22_HYD = %lf \n", S22_hyd);
				printf("S33_HYD = %lf \n", S33_hyd);
				printf("div_v = %lf \n", div_v);
				printf("c = %lf \n", c);
				printf("Le = %lf \n", Le);

			}

			S11 = stressVoigt->ve[0] + S11_hyd; S12 = stressVoigt->ve[2];
			S21 = stressVoigt->ve[3]; S22 = stressVoigt->ve[1]+S22_hyd;
			S33 = stressVoigt->ve[4] + S33_hyd;

			F11 = scni[i]->F_r->me[0][0]; F12 = scni[i]->F_r->me[0][1];
			F21 = scni[i]->F_r->me[1][0]; F22 = scni[i]->F_r->me[1][1];
			F33 = scni[i]->F_r->me[2][2];

			stressVoigt->ve[0] = S11*F11 + S12*F12;
			stressVoigt->ve[1] = S21*F21 + S22*F22; 
			stressVoigt->ve[2] = S11*F21 + S12*F22;
			stressVoigt->ve[3] = S21*F11 + S22*F12;
			stressVoigt->ve[4] = S33*F33;

		}else{
			fprintf(stderr,"dimension not yet supported");
		}


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
	M_FREE(invF);
	M_FREE(L);

}
	++call_count;

	return delta_t_min;
}
