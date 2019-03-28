
#include "internal_force_mooney.h"
#include "Deformation/velocity_grad.h"


static double epsilon_penalty = -1e6;

static double mu = 4000;
static double lambda = 1e5;
static double rho = 1000e-9;


//extern double * FINT_GLOBAL;
//extern double * RPEN_GLOBAL;




void internal_force_mooney(void *threadarg) {
	

	// Initial variables
	double P11,P12,P13,P21,P22,P23,P31,P32,P33;
	double Xp_n, Xp_n_1, Yp_n, Yp_n_1;
	internal_force_args * internal_force_struct = (internal_force_args*) threadarg;
	MAT * F_r = internal_force_struct->MP->F_n;
	MATERIAL_POINT * MP = internal_force_struct->MP;
	int num_neighbours = internal_force_struct->MP->num_neighbours;


	// Compute incremental deformation gradient
	defgrad_m(internal_force_struct->MP->inc_F, 
		internal_force_struct->MP->B, 
		internal_force_struct->MP->neighbours, 
		num_neighbours, internal_force_struct->inc_disp);


	// Compute total deformation gradient 
	m_mlt(internal_force_struct->MP->inc_F,
		internal_force_struct->MP->F_n,
		internal_force_struct->MP->stateNew->F);

	
	// Integration factor
	double intFactor = internal_force_struct->MP->volume*
	internal_force_struct->MP->INTEGRATION_FACTOR;

	// Material law 
	mooneyRivlin(internal_force_struct->MP->stressVoigt,
		internal_force_struct->MP->stateNew,
		internal_force_struct->materialParameters);




	// push forward piola kirchoff stress to Omega_n configuration
	sv_mlt(1.00/internal_force_struct->MP->Jn,
		internal_force_struct->MP->stressVoigt,
		internal_force_struct->MP->stressVoigt);


	P11 = internal_force_struct->MP->stressVoigt->ve[0];
	P22 = internal_force_struct->MP->stressVoigt->ve[1];
	P12 = internal_force_struct->MP->stressVoigt->ve[2];
	P21 = internal_force_struct->MP->stressVoigt->ve[3];


	internal_force_struct->MP->stressVoigt->ve[0] = P11*F_r->me[0][0] + P12*F_r->me[0][1];
	internal_force_struct->MP->stressVoigt->ve[1] = P21*F_r->me[1][0] + P22*F_r->me[1][1];
	internal_force_struct->MP->stressVoigt->ve[2] = P11*F_r->me[1][0] + P12*F_r->me[1][1];
	internal_force_struct->MP->stressVoigt->ve[3] = P21*F_r->me[0][0] + P22*F_r->me[0][1];



	// Get internal force for this material point
	__zero__(internal_force_struct->MP->fInt->ve,internal_force_struct->MP->fInt->max_dim);

	vm_mlt(internal_force_struct->MP->B,
		internal_force_struct->MP->stressVoigt,
		internal_force_struct->MP->fInt);


	// Assemble into global internal force vector
	for ( int k = 0 ; k < MP->num_neighbours; k++){

		int index = MP->neighbours->ive[k];
		internal_force_struct->FINT->ve[2*index] += intFactor * MP->fInt->ve[2*k];
		internal_force_struct->FINT->ve[2*index + 1] += intFactor *MP->fInt->ve[2*k+1];


	}


	// update_material point neighbours
	MP = update_material_point(MP, internal_force_struct->cells,
		internal_force_struct->XI_n_1, internal_force_struct->NODAL_MASS);


	// * PENALTY METHOD TO ENFORCE LINEAR DISPLACEMENT FIELD *// 

	// Find error in displacment field
	Xp_n_1 = MP->coords_n_1[0];
	Xp_n = MP->coords_n[0];
	Yp_n_1 = MP->coords_n_1[1];
	Yp_n = MP->coords_n[1];
	num_neighbours = MP->num_neighbours;
	MAT * XI_n_1 = internal_force_struct->XI_n_1;
	MAT * XI_n = internal_force_struct->XI_n;


	for ( int k = 0 ; k < num_neighbours; k++){

		int index = MP->neighbours->ive[k];

		// Find vectors dX_n and dX_n_1
		double delta_x_n[2] = {XI_n->me[index][0] - Xp_n,XI_n->me[index][1] - Yp_n};
		double delta_x_n_1[2] = {XI_n_1->me[index][0] - Xp_n_1,XI_n_1->me[index][1] - Yp_n_1};

		// Find error in displacement 
		double x_tilde[2] = {MP->inc_F->me[0][0]*(delta_x_n[0])
			+ MP->inc_F->me[0][1]*(delta_x_n[1]),
			MP->inc_F->me[1][0]*delta_x_n[0]
			+ MP->inc_F->me[1][1]*delta_x_n[1]};

		double norm_x = sqrt(pow(delta_x_n[0],2) + pow(delta_x_n[1],2));
		double e_x  = (delta_x_n_1[0] - x_tilde[0])/norm_x;
		double e_y =  (delta_x_n_1[1] - x_tilde[1])/norm_x ;

		// assemble forces 
		internal_force_struct->RPEN->ve[2*index] += epsilon_penalty*MP->shape_function->phi->ve[k]*e_x*intFactor;
		internal_force_struct->RPEN->ve[2*index+1] += epsilon_penalty*MP->shape_function->phi->ve[k]*e_y*intFactor;

		}

		return ;
	}


//{
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
	// double dt = internal_force_struct->dt;

	// __zero__(FINT->ve,FINT->max_dim);
	// __zero__(NODAL_MASS->ve,NODAL_MASS->max_dim);
	// __zero__(RPEN->ve,RPEN->max_dim);


	// for ( int k = 0 ; k < mat_points->max_dim; k++)
	// {


	// 	i = mat_points->ive[k];
	// 	F = material_points->MP[i]->stateNew->F;
	// 	Fdot = material_points->MP[i]->stateNew->Fdot;
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


	// 	// // Find Fdot
	// 	// m_sub(material_points->MP[i]->stateNew->F,material_points->MP[i]->stateOld->F,
	// 	// 	material_points->MP[i]->stateNew->m_temp1);
	// 	// sm_mlt(1.00/dt,material_points->MP[i]->stateNew->m_temp1,material_points->MP[i]->stateNew->Fdot);
	// 	// Find L
	// 	// m_mlt(stateNew->invF,stateNew->Fdot,stateNew->L);

	// 	// velocity_grad(stateNew, stateOld, dt,0.5);
	// 	// m_inverse_small(stateNew->F,stateNew->invF);
	// 	// stateNew->Jacobian = determinant(stateNew->F);


	// 	// double div_v = stateNew->div_v;


	// 	// double b1 = 0.06;
	// 	// double b2 = 1.44;

	// 	// double Le = sqrt(material_points->MP[i]->volume);
	// 	// //Le = 2;
	// 	// double c = sqrt(((lambda+2*mu)/rho));
	// 	// double P_b1 = 0;

	// 	// P_b1 = b1*div_v*rho*Le*c;
	// 	// double P_b2 = 0;
	// 	// if ( div_v < 0 ){
	// 	// 	P_b2 = Le*rho*(b2*b2*Le*div_v*div_v) ;
	// 	// }

	// 	// double P11_hydr = (stateNew->Jacobian)*(P_b1 + P_b2)*stateNew->invF->me[0][0];
	// 	// double P22_hydr = (stateNew->Jacobian)*(P_b1 + P_b2)*stateNew->invF->me[1][1];
	// 	// // Compute 


	// 	double P11_hydr = 0;
	// 	double P22_hydr = 0;

	// 	// Integration factor
	// 	double intFactor = material_points->MP[i]->volume*
	// 	material_points->MP[i]->INTEGRATION_FACTOR;


	// 	// Material law 
	// 	mooneyRivlin(stressVoigt,material_points->MP[i]->stateNew,materialParameters);


	// 	__zero__(fInt->ve,fInt->max_dim);


	// 	// push forward piola kirchoff stress to Omega_n configuration
	// 	sv_mlt(1.00/material_points->MP[i]->Jn,stressVoigt,stressVoigt);

	// 	P11 = stressVoigt->ve[0]+P11_hydr;
	// 	P22 = stressVoigt->ve[1]+P22_hydr;
	// 	P12 = stressVoigt->ve[2];
	// 	P21 = stressVoigt->ve[3];

	// 	stressVoigt->ve[0] = P11*F_r->me[0][0] + P12*F_r->me[0][1];
	// 	stressVoigt->ve[1] = P21*F_r->me[1][0] + P22*F_r->me[1][1];
	// 	stressVoigt->ve[2] = P11*F_r->me[1][0] + P12*F_r->me[1][1];
	// 	stressVoigt->ve[3] = P21*F_r->me[0][0] + P22*F_r->me[0][1];


	// 	vm_mlt(B,stressVoigt,fInt);



	// 	// Assemble internal force vector
	// 	for ( int k = 0 ; k < num_neighbours; k++){

	// 		int index = neighbours->ive[k];

	// 		FINT->ve[2*index] += intFactor * fInt->ve[2*k];
	// 		FINT->ve[2*index + 1] += intFactor *fInt->ve[2*k+1];


	// 	}


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
	// 		RPEN->ve[2*index] += epsilon_penalty*material_points->MP[i]->shape_function->phi->ve[k]*e_x*intFactor;
	// 		RPEN->ve[2*index+1] += epsilon_penalty*material_points->MP[i]->shape_function->phi->ve[k]*e_y*intFactor;

	// 		}




	// 	} //end of loop over points




	// 	return ;
	// }


