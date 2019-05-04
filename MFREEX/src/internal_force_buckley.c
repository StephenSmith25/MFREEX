
#include "internal_force_buckley.h"



static double epsilon_penalty = 50;  // NORMAL VLAUE = -50
static double xi =0.01;

#define HOURGLASS_CONTROL 
//#define VISCOUS_HOURGLASS
#define STIFFNESS_HOURGLASS

void internal_force_buckley(void *threadarg){



	// Initial variables
	double P11,P12,P13,P21,P22,P23,P31,P32,P33;
	double Xp_n, Xp_n_1, Yp_n, Yp_n_1;
	internal_force_args * internal_force_struct = (internal_force_args*) threadarg;
	MAT * F_r = internal_force_struct->MP->F_n;
	MATERIAL_POINT * MP = internal_force_struct->MP;
	int num_neighbours = internal_force_struct->MP->num_neighbours;
	double DT = internal_force_struct->dt;

	// Compute incremental deformation gradient
	defgrad_m(internal_force_struct->MP->inc_F, 
		internal_force_struct->MP->B, 
		internal_force_struct->MP->neighbours, 
		num_neighbours, internal_force_struct->inc_disp);


	// Compute total deformation gradient 
	m_mlt(internal_force_struct->MP->inc_F,
		internal_force_struct->MP->F_n,
		internal_force_struct->MP->stateNew->F);


	/* ------------------------------------------*/
	/* -----------------Deformation--------------*/
	/* ------------------------------------------*/
	/*  Find Rate of deformation tensors at n+h*/
	velocity_grad(MP->stateNew,MP->stateOld, DT, 0.500);

	// Find inverse deformation gradient at n+1
	m_inverse_small(MP->stateNew->F, MP->stateNew->invF);


	// Find Jacobian at n+1
	MP->stateNew->Jacobian = determinant(MP->stateNew->F);

	//------------------------------------------//
	//        Update Polar Decomposition        //
	//------------------------------------------//
	update_Polar_Decomposition(MP->stateNew, MP->stateOld, DT);

	// remove rotation from D to get d
	un_rotate_tensor(MP->stateNew->D, MP->stateNew->R, 
		MP->stateNew->m_temp1, MP->stateNew->d);

	/* ------------------------------------------*/
	/* -----------------Damping -----------------*/
	/* ------------------------------------------*/
	double rho = 1380e-9;
	double b1 = 0;
	double b2 = 1.44;
	double Le =1;
	double Cd = 1400;
	double div_v = MP->stateNew->div_v;
	double qv =  rho*Le*b1*Cd * div_v;
	if ( div_v < 0)
	{
		qv += rho*Le*(b2 * (Le/1000) * pow(div_v,2)) ;
	}
	qv = 0;

	// Integration factor
	double intFactor = MP->volume*MP->INTEGRATION_FACTOR;


	// -------------------------------------------------//
	//          Get Voigt Piola Kirchoff stress         //
	// -------------------------------------------------//

	//Material law 
	buckleyStress(MP->stateNew,
		MP->stateOld,
		internal_force_struct->materialParameters,
		DT);


	internal_force_struct->sigma->ve[0] = (MP->stateNew->sigma->me[0][0])/1e6 + qv;
	internal_force_struct->sigma->ve[1] = (MP->stateNew->sigma->me[1][1])/1e6 + qv;
	internal_force_struct->sigma->ve[2] = MP->stateNew->sigma->me[0][1]/1e6;
	internal_force_struct->sigma->ve[3] = (MP->stateNew->sigma->me[2][2])/1e6 + qv;



	/*  Internal force vectors */
	gMat(internal_force_struct->G,MP->stateNew->invF,1);
	mv_mlt(internal_force_struct->G,internal_force_struct->sigma,MP->stressVoigt);
	sv_mlt(MP->stateNew->Jacobian,MP->stressVoigt,MP->stressVoigt);


	// push forward piola kirchoff stress to Omega_n configuration
	sv_mlt(1.00/MP->Jn,MP->stressVoigt,MP->stressVoigt);

	P11 = MP->stressVoigt->ve[0];
	P22 = MP->stressVoigt->ve[1];
	P12 = MP->stressVoigt->ve[2];
	P21 = MP->stressVoigt->ve[3];
	P33 = MP->stressVoigt->ve[4];

	MP->stressVoigt->ve[0] = P11*F_r->me[0][0] + P12*F_r->me[0][1];
	MP->stressVoigt->ve[1] = P21*F_r->me[1][0] + P22*F_r->me[1][1];
	MP->stressVoigt->ve[2] = P11*F_r->me[1][0] + P12*F_r->me[1][1];
	MP->stressVoigt->ve[3] = P21*F_r->me[0][0] + P22*F_r->me[0][1];
	MP->stressVoigt->ve[4] = P33*F_r->me[2][2];



	// -------------------------------------------------//
	//           	 Internal force vector               //
	// -------------------------------------------------//

	// Get internal force for this material point
	__zero__(internal_force_struct->MP->fInt->ve,internal_force_struct->MP->fInt->max_dim);

	vm_mlt(MP->B,MP->stressVoigt,MP->fInt);



	// Assemble into global internal force vector
	for ( int k = 0 ; k < MP->num_neighbours; k++){

		int index = MP->neighbours->ive[k];
		internal_force_struct->FINT->ve[2*index] += intFactor * MP->fInt->ve[2*k];
		internal_force_struct->FINT->ve[2*index + 1] += intFactor *MP->fInt->ve[2*k+1];


	}


	// -------------------------------------------------//
	//                  State storage                   //
	// -------------------------------------------------//

	// n+1 --->>>> n 
	m_copy(MP->stateNew->F,MP->stateOld->F);
	m_copy(MP->stateNew->V,MP->stateOld->V);
	m_copy(MP->stateNew->R,MP->stateOld->R);
	m_copy(MP->stateNew->invF,MP->stateOld->invF);			
	MP->stateOld->Jacobian = MP->stateNew->Jacobian;
	MP->stateOld->div_v = MP->stateNew->div_v;



	// -------------------------------------------------//
	//             Penalty stabilisation                 //
	// -------------------------------------------------//

	// * PENALTY METHOD TO ENFORCE LINEAR DISPLACEMENT FIELD *// 

	MAT * XI_n_1 = internal_force_struct->XI_n_1;
	MAT * XI_n = internal_force_struct->XI_n;
	/* UPDATE MATERIAL POINT COORDINATES*/
	for ( int i = 0 ; i < MP->num_neighbours ; i++)
	{	
		int index = MP->neighbours->ive[i];

		for ( int k = 0 ; k < 2 ; k++)
		{
			if ( i == 0)
			{
				MP->coords_n_1[k] = 0;	
			}
			MP->coords_n_1[k] += MP->shape_function->phi->ve[i] * XI_n_1->me[index][k];
		}
	}

#ifdef HOURGLASS_CONTROL


	VEC * velocity = internal_force_struct->velocity;
	// Find error in displacment field
	Xp_n_1 = MP->coords_n_1[0];
	Xp_n = MP->coords_n[0];
	Yp_n_1 = MP->coords_n_1[1];
	Yp_n = MP->coords_n[1];
	num_neighbours = MP->num_neighbours;
	VEC * mass = internal_force_struct->mass;



#ifdef VISCOUS_HOURGLASS
	double v_p_x = 0;
	double v_p_y = 0;
	/* UPDATE MATERIAL POINT COORDINATES*/
	for ( int i = 0 ; i < MP->num_neighbours ; i++)
	{	
		int index = MP->neighbours->ive[i];

		v_p_x += MP->shape_function->phi->ve[i] * velocity->ve[2*index];
		v_p_y += MP->shape_function->phi->ve[i] * velocity->ve[2*index+1];

	}

#endif


	for ( int k = 0 ; k < num_neighbours; k++){


		double phi = MP->shape_function->phi->ve[k];
		int index = MP->neighbours->ive[k];
		// Find vectors dX_n and dX_n_1
		double X_j_p[2] = {XI_n->me[index][0] - Xp_n,XI_n->me[index][1] - Yp_n};
		double x_j_p[2] = {XI_n_1->me[index][0] - Xp_n_1,XI_n_1->me[index][1] - Yp_n_1};
		// Find error in displacement 
		double x_tilde[2] = {MP->inc_F->me[0][0]*(X_j_p[0])
			+ MP->inc_F->me[0][1]*(X_j_p[1]),
			MP->inc_F->me[1][0]*X_j_p[0]
			+ MP->inc_F->me[1][1]*X_j_p[1]};

		double norm_X = sqrt(pow(X_j_p[0],2) + pow(X_j_p[1],2));
		double norm_x = sqrt(pow(x_j_p[0],2) + pow(x_j_p[1],2));

		double epsilon[2] = {x_tilde[0] -  x_j_p[0] ,x_tilde[1] -  x_j_p[1] };


#ifdef STIFFNESS_HOURGLASS


		// ******************************************** //
		//			Stiffness hour glass control	   //
		// ********************************************//

		double e_x  = epsilon[0] /norm_X;
		double e_y = epsilon[1] /norm_X;

		// assemble forces 
		internal_force_struct->RPEN->ve[2*index] += epsilon_penalty*MP->shape_function->phi->ve[k]*e_x*intFactor;
		internal_force_struct->RPEN->ve[2*index+1] += epsilon_penalty*MP->shape_function->phi->ve[k]*e_y*intFactor;

#endif
		// ******************************************** //
		//			 Viscous hour glass control	       //
		// ********************************************//
#ifdef VISCOUS_HOURGLASS

		double mass = internal_force_struct->mass->ve[index];

		double h = MP->r_cutoff;

		// Find error based on deformation gradient
		double error_mag = sqrt(pow((x_tilde[0] -  x_j_p[0] ),2) + pow(x_tilde[1] -  x_j_p[1],2));
		double error_mag_normalised = error_mag/norm_X;
		double unit_vectors_x[2] = {x_j_p[0]/norm_x, x_j_p[1]/norm_x};

		// epsilon / || X || *  x / || x || 
		double viscous_ep_x = unit_vectors_x[0] * error_mag_normalised;
		double viscous_ep_y = unit_vectors_x[1] * error_mag_normalised;


		// check if need to apply viscous control
		double v_j_x = velocity->ve[2*index];
		double v_j_y = velocity->ve[2*index+1];

		// v_p - v_j
		double v_rel_x = v_p_x - v_j_x;
		double v_rel_y = v_p_y - v_j_y;
		// v_pj dot dx 
		double v_dot_x = v_rel_x*x_j_p[0] + v_rel_y*x_j_p[1]  ;

		// Epsilon_pj dot x_pj
		double ep_dot_x = epsilon[0]*x_j_p[0] + epsilon[1]*x_j_p[1];

		double flag = v_dot_x * ep_dot_x;

		if ( flag <= 0)
		{
				internal_force_struct->RPEN->ve[2*index] += mass*xi*v_dot_x*phi*viscous_ep_x*h*intFactor*(Cd*1000);
				internal_force_struct->RPEN->ve[2*index+1] += mass*xi*v_dot_x*phi*viscous_ep_y*h*intFactor*(Cd*1000);

		}

	
#endif
	

		}
#endif



		return ;
	}


	// // Internal forces
	// double P11,P12,P13,P21,P22,P23,P31,P32,P33;

	// MAT * F ;
	// MAT * Fdot;
	// VEC * stressVoigt ;
	// MAT * B ;
	// IVEC * neighbours;
	// MAT * F_r;
	// VEC * fInt;
	// double Xp_n, Xp_n_1, Yp_n, Yp_n_1;
	// int i = 0;

	// internal_force_args * internal_force_struct = (internal_force_args*) threadarg;

	// // Get inputs from struct 
	// VEC * FINT = internal_force_struct->FINT;
	// VEC * NODAL_MASS = internal_force_struct->NODAL_MASS;
	// VEC * RPEN = internal_force_struct->RPEN;
	// IVEC * mat_points = internal_force_struct->mat_points;
	// VEC * materialParameters = internal_force_struct->materialParameters;
	// VEC * inc_disp = internal_force_struct->inc_disp;
	// CELLS * cells = internal_force_struct->cells;
	// MATERIAL_POINTS * material_points = internal_force_struct->material_points;
	// MAT * XI_n = internal_force_struct->XI_n;
	// MAT * XI_n_1 = internal_force_struct->XI_n_1;
	// double DT = internal_force_struct->dt;
	// MAT * G = internal_force_struct->G;
	// VEC * sigma = internal_force_struct->sigma;

	// __zero__(FINT->ve,FINT->max_dim);
	// __zero__(NODAL_MASS->ve,NODAL_MASS->max_dim);
	// __zero__(RPEN->ve,RPEN->max_dim);


	// double epsilon_penalty = 0;




	// for ( int k = 0 ; k < mat_points->max_dim; k++)
	// {



	// 	int i = mat_points->ive[k];

	// 	state_variables *stateNew = material_points->MP[i]->stateNew;
	// 	state_variables *stateOld = material_points->MP[i]->stateOld;


	// 	F = material_points->MP[i]->stateNew->F;
	// 	B = material_points->MP[i]->B;
	// 	neighbours = material_points->MP[i]->neighbours;
	// 	F_r = material_points->MP[i]->F_n;
	// 	stressVoigt = material_points->MP[i]->stressVoigt;
	// 	fInt = material_points->MP[i]->fInt;
	// 	int num_neighbours = material_points->MP[i]->num_neighbours;




	// 	// Compute incremental deformation gradient
	// 	defgrad_m(material_points->MP[i]->inc_F, B, neighbours, num_neighbours, inc_disp);

	// 	// compute total deformation gradient 
	// 	m_mlt(material_points->MP[i]->inc_F,
	// 		material_points->MP[i]->F_n,material_points->MP[i]->stateNew->F);


	// 	/* ------------------------------------------*/
	// 	/* -----------------Deformation--------------*/
	// 	/* ------------------------------------------*/
	// 	/*  Find Rate of deformation tensors at n+h*/
	// 	velocity_grad(stateNew,stateOld, DT, 0.500);

	// 	// Find inverse deformation gradient at n+1
	// 	m_inverse_small(stateNew->F, stateNew->invF);


	// 	// Find Jacobian at n+1
	// 	stateNew->Jacobian = determinant(stateNew->F);


	// 	//------------------------------------------//
	// 	//        Update Polar Decomposition        //
	// 	//------------------------------------------//
	// 	update_Polar_Decomposition(stateNew, stateOld, DT);

	// 	// remove rotation from D to get d
	// 	un_rotate_tensor(stateNew->D, stateNew->R, 
	// 			stateNew->m_temp1, stateNew->d);

	// 	/* ------------------------------------------*/
	// 	/* -----------------Damping -----------------*/
	// 	/* ------------------------------------------*/
	// 	double rho = 1380e-9;

	// 	double b1 = 0.06;
	// 	double b2 = 1.44;
	// 	double Le =0;
	// 	double Cd = 1400;
	// 	double div_v = stateNew->div_v;
	// 	double qv =  rho*Le*b1*Cd * div_v;
	// 	if ( div_v < 0)
	// 	{
	// 		qv += rho*Le*(b2 * (Le/1000) * pow(div_v,2)) ;
	// 	}


	// 	// Integration factor
	// 	double intFactor = material_points->MP[i]->volume*
	// 	material_points->MP[i]->INTEGRATION_FACTOR;


	// 	// Material law 
	// 	buckleyStress(stateNew,stateOld,materialParameters,DT);
	// 	//m_foutput(stdout, stateNew->sigma);


	// 	__zero__(fInt->ve,fInt->max_dim);



	// 	//-------------------------------------------------//
	// 	//         Get Voigt Piola Kirchoff stress         //
	// 	//-------------------------------------------------//
	// 	// v_zero(sigma);
	// 	// v_zero(stressVoigt);
	// 	// m_zero(G);

	// 	sigma->ve[0] = (stateNew->sigma->me[0][0])/1e6 + qv;
	// 	sigma->ve[1] = (stateNew->sigma->me[1][1])/1e6 + qv;
	// 	sigma->ve[2] = stateNew->sigma->me[0][1]/1e6;
	// 	sigma->ve[3] = (stateNew->sigma->me[2][2])/1e6 + qv;



	// 	/*  Internal force vectors */
	// 	gMat(G,stateNew->invF,1);
	// 	mv_mlt(G,sigma,stressVoigt);
	// 	sv_mlt(stateNew->Jacobian,stressVoigt,stressVoigt);


	// 	// push forward piola kirchoff stress to Omega_n configuration
	// 	sv_mlt(1.00/material_points->MP[i]->Jn,stressVoigt,stressVoigt);

	// 	P11 = stressVoigt->ve[0];
	// 	P22 = stressVoigt->ve[1];
	// 	P12 = stressVoigt->ve[2];
	// 	P21 = stressVoigt->ve[3];
	// 	P33 = stressVoigt->ve[4];

	// 	stressVoigt->ve[0] = P11*F_r->me[0][0] + P12*F_r->me[0][1];
	// 	stressVoigt->ve[1] = P21*F_r->me[1][0] + P22*F_r->me[1][1];
	// 	stressVoigt->ve[2] = P11*F_r->me[1][0] + P12*F_r->me[1][1];
	// 	stressVoigt->ve[3] = P21*F_r->me[0][0] + P22*F_r->me[0][1];
	// 	stressVoigt->ve[4] = P33*F_r->me[2][2];

	// 	vm_mlt(B,stressVoigt,fInt);

	// 	// Assemble internal force vector
	// 	for ( int l = 0 ; l < num_neighbours; l++){

	// 		int index = neighbours->ive[l];

	// 		FINT->ve[2*index] += intFactor * fInt->ve[2*l];
	// 		FINT->ve[2*index + 1] += intFactor *fInt->ve[2*l+1];


	// 	}

	// 	// n+1 --->>>> n 
	// 	m_copy(stateNew->F,stateOld->F);
	// 	m_copy(stateNew->V,stateOld->V);
	// 	m_copy(stateNew->R,stateOld->R);
	// 	m_copy(stateNew->invF,stateOld->invF);			
	// 	stateOld->Jacobian = stateNew->Jacobian
	// 	stateOld->div_v = stateNew->div_v;



	// 	// update_material point neighbours

	// 	material_points->MP[i] = update_material_point(material_points->MP[i], cells,
	// 		 XI_n_1, NODAL_MASS);

	// 	// Find error in displacment field
	// 	Xp_n_1 = material_points->MP[i]->coords_n_1[0];
	// 	Xp_n = material_points->MP[i]->coords_n[0];
	// 	Yp_n_1 = material_points->MP[i]->coords_n_1[1];
	// 	Yp_n = material_points->MP[i]->coords_n[1];


	// 	num_neighbours = material_points->MP[i]->num_neighbours;
	// 	neighbours = material_points->MP[i]->neighbours;


	// 	for ( int k = 0 ; k < num_neighbours; k++){

	// 		int index = neighbours->ive[k];

	// 		// Find vectors dX_n and dX_n_1
	// 		double delta_x_n[2] = {XI_n->me[index][0] - Xp_n,XI_n->me[index][1] - Yp_n};
	// 		double delta_x_n_1[2] = {XI_n_1->me[index][0] - Xp_n_1,XI_n_1->me[index][1] - Yp_n_1};

	// 		// Find error in displacement 
	// 		double x_tilde[2] = {material_points->MP[i]->inc_F->me[0][0]*(delta_x_n[0])
	// 				+ material_points->MP[i]->inc_F->me[0][1]*(delta_x_n[1]),
	// 				material_points->MP[i]->inc_F->me[1][0]*delta_x_n[0]
	// 				+ material_points->MP[i]->inc_F->me[1][1]*delta_x_n[1]};

	// 		double norm_x = sqrt(pow(delta_x_n[0],2) + pow(delta_x_n[1],2));
	// 		double e_x  = (delta_x_n_1[0] - x_tilde[0])/norm_x;
	// 		double e_y =  (delta_x_n_1[1] - x_tilde[1])/norm_x ;

	// 		// assemble forces 
	// 		RPEN->ve[2*index] += epsilon_penalty*material_points->MP[i]->shape_function->phi->ve[k]*e_x;
	// 		RPEN->ve[2*index+1] += epsilon_penalty*material_points->MP[i]->shape_function->phi->ve[k]*e_y;

	// 		}




	// 	} //end of loop over points




	// 	return ;