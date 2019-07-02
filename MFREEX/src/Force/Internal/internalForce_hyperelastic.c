#include "Force/Internal/internalForce_hyperelastic.h"
static int call_count ;
static int print_count;


double internalForce_hyperelastic(VEC * Fint, MATERIAL_POINT * MP, VEC * disp, VEC * velocity,
 	VEC * matParams, int (*mat_func_ptr)(VEC *, MAT*,VEC*), double t_n_1){


		// Material point 
		state_variables * stateNew = MP->stateNew;


		// time step calculation
		// 
		double delta_t_min = 1000;


		MAT * F = stateNew->F;
		 MAT * Fdot = stateNew->Fdot;
		// MAT * L = stateNew->L;
		 VEC * stressVoigt = MP->stressVoigt;
		// MAT * invF = stateNew->invF;
		// int dim = F->m;



		// // stress tensor
		// double S11,S12,S13,S21,S22,S23,S31,S32,S33;
		// double F11,F12,F13,F21,F22,F23,F31,F32,F33;
		// double invF11,invF12,invF13,invF21,invF22,invF23,invF31,invF32,invF33;
		// double Co11,Co12,Co13,Co21,Co22,Co23,Co31,Co32,Co33;
		// double L11,L12,L13,L21,L22,L23,L31,L32,L33;
		// double div_v;
		// double delta_t_min_i = 1000;
		MAT * B ;
		IVEC * neighbours;
		MAT * F_r;
		VEC * fInt;

		int num_neighbours = MP->num_neighbours;


		/*  Find deformation gradient */
		B = MP->B;
		neighbours = MP->neighbours;
		F_r = MP->F_n;
		fInt = MP->fInt;


		/*  Find deformation gradient */
		//get_defgrad(F, B, neighbours,F_r, disp);

		defgrad_m(MP->inc_F, 
		MP->B, 
		MP->neighbours, 
		num_neighbours, 
		disp);



		// Compute total deformation gradient 
		m_mlt(MP->inc_F,
		MP->F_n,
		F);

	




		get_dot_defgrad(Fdot,B, neighbours,F_r,velocity);

	

		

		//double Jacobian ;


		// if ( dim == 2)
		// {	

		// 	F11 = F->me[0][0]; F12 = F->me[0][1];
		// 	F21 = F->me[1][0]; F22 = F->me[1][1];


		// 	Jacobian = F11*F22 - F12*F21;

		// 	Co11 = F22;
		// 	Co12 = -F21;
		// 	Co21 = -F12;
		// 	Co22 = F11;
		// 	invF11 = Co11/Jacobian; invF12 = Co21/Jacobian;
		// 	invF21 = Co12/Jacobian; invF22 = Co22/Jacobian;


		// 	// Find velocity gradient

		// 	L11 = Fdot->me[0][0]*invF11 + Fdot->me[0][1]*invF21; L12 = Fdot->me[0][0]*invF12 + Fdot->me[0][1]*invF22;
		// 	L21 = Fdot->me[1][0]*invF11 + Fdot->me[1][1]*invF21; L22 = Fdot->me[1][0]*invF12 + Fdot->me[1][1]*invF22;



		// 	div_v = L11 + L22;
		// }else if ( dim == 3)
		// {

		// 	m_inverse_small(F,invF);

		// 	m_mlt(Fdot,invF,L);

		// 	Jacobian = determinant(F);

		// 	// Find velocity gradient
		// 	div_v = L->me[0][0] + L->me[1][1] + L->me[2][2];	
		// }


		/* Integration parameter */
		double intFactor = MP->volume;



		// //Find 1D frequency bounds

		// double B11 = 0;
		// double B22 = 0;
		// double B33 = 0;

		// double fro_b = 0;

		// double MaxB = -1;
		// double rho =  1800e-9;
		// double num_neighbours_J = 0;
		// double m_j;
		// double volume_I = intFactor;
		// for ( int k = 0 ; k < (num_neighbours) ; k++)
		// {	
		// 	int index_J = neighbours->ive[k];
		// 	num_neighbours_J = MPS[index_J]->neighbours->max_dim;
		// 	volume_I = MP->volume;
		// 	m_j = rho * MPS[index_J]->volume;


		// 	fro_b += (volume_I*num_neighbours_J/m_j)*B->me[0][2*k]*B->me[0][2*k];
		// 	fro_b += (volume_I*num_neighbours_J/m_j)*B->me[1][2*k+1]*B->me[1][2*k+1];
		// 	fro_b += (volume_I*num_neighbours_J/m_j)*B->me[2][2*k]*B->me[2][2*k];
		// 	fro_b += (volume_I*num_neighbours_J/m_j)*B->me[3][2*k+1]*B->me[3][2*k+1];

		// 	// if ( is_axi == 1){
		// 	// 	m_j = rho * scni[index_J]->area*2*PI*scni[index_J]->r;
		// 	// }

		// 	//B11 += (volume_I*num_neighbours_J/m_j)*(B->me[0][2*k]*B->me[0][2*k]);
		// 	//B22 += (volume_I*num_neighbours_J/m_j)*(B->me[1][2*k+1]*B->me[1][2*k+1]);

		// 	if ( is_axi == 1)
		// 	{
		// 		//B33 += (volume_I*num_neighbours_J/m_j)*(B->me[4][2*k]*B->me[4][2*k]);
		// 	fro_b += (volume_I*num_neighbours_J/m_j)*B->me[4][2*k]*B->me[4][2*k];

		// 	}

		// }


		// //MaxB = max(B11,B22);

		// MaxB = fro_b;

		// double omega_max = sqrt((lambda + 2*mu)*fro_b);



		// // if ( is_axi == 1)
		// // {
		// // 	MaxB = max(MaxB,B33);

		// // }

		// double b1 = 0;
		// double b2 = 0;

		// double Le = sqrt(MP->volume);
		// //Le = 2;
		// double c = sqrt(((lambda+2*mu)/rho));
		// double P_b1 = 0;
		// //P_b1 = b1*div_v*rho*Le*c;
		// double eta = b1;
		// double P_b2 = 0;
		// if ( div_v < 0 ){
		// 	P_b2 = Le*rho*(b2*b2*Le*div_v*div_v) -  b1*div_v*rho*Le*c;
		// 	eta -= b2*b2*Le*(1/c)*div_v;
		// }

		// double delta_t = (2.00/omega_max)*(sqrt(1+eta*eta)-eta); 
		// //double delta_t = (2.00/sqrt(((lambda + 2*mu)*MaxB)))*(sqrt(1+eta*eta)-eta);
		// if ( delta_t <delta_t_min_i )
		// {
		// 	delta_t_min_i = delta_t;
		// }


		// Get stress 
		mat_func_ptr(stressVoigt,F,matParams);


		__zero__(fInt->ve,fInt->max_dim);


		// // push forward stress to new reference configuration
		// if ( dim_v == 4)
		// {

		// 	// double S11_hyd = (Jacobian)*(P_b1+P_b2)*invF11;
		// 	// double S22_hyd = (Jacobian)*(P_b1+P_b2)*invF22;


		// 	// S11 = stressVoigt->ve[0]+S11_hyd; S12 = stressVoigt->ve[2];
		// 	// S21 = stressVoigt->ve[3]; S22 = stressVoigt->ve[1]+S22_hyd;

		// 	// F11 = scni[i]->F_r->me[0][0]; F12 = scni[i]->F_r->me[0][1];
		// 	// F21 = scni[i]->F_r->me[1][0]; F22 = scni[i]->F_r->me[1][1];

		// 	// stressVoigt->ve[0] = S11*F11 + S12*F12;
		// 	// stressVoigt->ve[1] = S21*F21 + S22*F22; 
		// 	// stressVoigt->ve[2] = S11*F21 + S12*F22;
		// 	// stressVoigt->ve[3] = S21*F11 + S22*F12;

		// }else if ( dim_v == 5)
		// {
		// 	// double S11_hyd = (Jacobian)*(P_b1+P_b2)*invF->me[0][0];
		// 	// double S22_hyd = (Jacobian)*(P_b1+P_b2)*invF->me[1][1];
		// 	// double S33_hyd = (Jacobian)*(P_b1+P_b2)*invF->me[2][2];

		// 	// if (( i == 0) && (call_count % 1000 == 0))
		// 	// {
		// 	// 	printf("S11_HYD = %lf \n", S11_hyd);
		// 	// 	printf("S22_HYD = %lf \n", S22_hyd);
		// 	// 	printf("S33_HYD = %lf \n", S33_hyd);
		// 	// 	printf("div_v = %lf \n", div_v);
		// 	// 	printf("c = %lf \n", c);
		// 	// 	printf("Le = %lf \n", Le);

		// 	// }

		// 	// S11 = stressVoigt->ve[0] + S11_hyd; S12 = stressVoigt->ve[2];
		// 	// S21 = stressVoigt->ve[3]; S22 = stressVoigt->ve[1]+S22_hyd;
		// 	// S33 = stressVoigt->ve[4] + S33_hyd;

		// 	// F11 = scni[i]->F_r->me[0][0]; F12 = scni[i]->F_r->me[0][1];
		// 	// F21 = scni[i]->F_r->me[1][0]; F22 = scni[i]->F_r->me[1][1];
		// 	// F33 = scni[i]->F_r->me[2][2];


		// 	// stressVoigt->ve[0] = S11*F11 + S12*F12;
		// 	// stressVoigt->ve[1] = S21*F21 + S22*F22; 
		// 	// stressVoigt->ve[2] = S11*F21 + S12*F22;
		// 	// stressVoigt->ve[3] = S21*F11 + S22*F12;
		// 	// stressVoigt->ve[4] = S33*F33;

		// }else{
		// 	fprintf(stderr,"dimension not yet supported");
		// }


		vm_mlt(B,stressVoigt,fInt);

		for ( int k = 0 ; k < num_neighbours; k++){
			Fint->ve[2*neighbours->ive[k]] += intFactor * fInt->ve[2*k];
			Fint->ve[2*neighbours->ive[k]+1] += intFactor *fInt->ve[2*k+1];
		}



	return delta_t_min;
}
