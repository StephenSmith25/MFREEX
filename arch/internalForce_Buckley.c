#include "Force/Internal/internalForce_Buckley.h"

static int call_count;
static int print_count = 1;
double internalForce_Buckley(VEC * Fint, SCNI_OBJ * scni_obj, VEC * disp, VEC * velocity,
	VEC * matParams,VEC * critLambdaParams, state_Buckley ** stateNew, state_Buckley ** stateOld,
	 int is_axi, int dim, double deltat, double t_n_1)
{


	__zero__(Fint->ve,Fint->max_dim);

	double Kb = matParams->ve[9];
	double Gb = matParams->ve[10];


	double mu = Gb;
	double lambda = Kb - (2.00/3.00)*mu;
	int dim_piola = 0;
	int dim_strain = 0;
	int dim_cauchy = 0;


	// check if problem is axisymmetric
	if ( is_axi == 1){
		dim_piola = 5;
		dim_strain = 3;
		dim_cauchy = 4;
	}else{
		dim_piola = dim*dim;
		dim_strain = dim;
		dim_cauchy = (dim*dim) - (dim -1);
	}

	// loop over all integration points

	SCNI ** scni = scni_obj->scni;
	int num_int_points = scni_obj->num_points;

	double averaging = 0;
	// time step calculation
	double delta_t_min = 1000;

	char  filename[50];


	// set number of threads
	omp_set_num_threads(8);

	int i;
#pragma omp parallel 
{


		// Internal force varaibles
		VEC * stressVoigt = v_get(dim_piola);
		VEC * sigma = v_get(dim_cauchy);
		VEC * fIntTemp = v_get(Fint->max_dim);
		double S11,S12,S13,S21,S22,S23,S31,S32,S33;
		double F11,F12,F13,F21,F22,F23,F31,F32,F33;
		// Gmat matrix
		MAT * G = m_get(dim_piola,dim_cauchy);

		// Time step
		double delta_t_min_i = 1000;

	
		MAT * B;
		IVEC * neighbours;
		MAT * F_r;
		MAT * Sc_n_1;
		MAT * Sb_n_1;


#pragma omp for nowait schedule(dynamic,2) 
		for(i = 0 ; i < num_int_points ; i++){



			/* ------------------------------------------*/
			/* -----------------Deformation--------------*/
			/* ------------------------------------------*/
		

			//------------------------//
			//        Def Grad        //
			//------------------------//

			/*  Find deformation gradient */
			B = scni[i]->B;
			neighbours = scni[i]->sfIndex;
			F_r = scni[i]->F_r;

			int num_neighbours = neighbours->max_dim;
			
			/*  Find deformation gradient */
			get_defgrad(stateNew[i]->F, B, neighbours,F_r,disp);

			buckleyStress(stateNew[i],stateOld[i],matParams,critLambdaParams,deltat,i,is_axi);

			/* ------------------------------------------*/
			/* --------------New time increment----------*/
			/* ------------------------------------------*/


			// Find 1D frequency bounds
			//time step calculation
			double rho = 1380e-9;
			double B11 = 0;
			double B22 = 0;
			double B33 = 0;
			double MaxB = -1;
			double num_neighbours_J = 0;
			double m_j;
			double volume_I = scni[i]->area;
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
				if ( is_axi == 1)
				{
					MaxB = max(MaxB,B33);
				}

			/* ------------------------------------------*/
			/* ---------------Bulk Damping---------------*/
			/* ------------------------------------------*/

			double div_v = stateNew[i]->div_v;
			double Jacobian = stateNew[i]->Jacobian;



			double lambda = stateNew[i]->lambda/1e6;
			double mu = stateNew[i]->mu/1e6;
			double b1 = 0;
			double b2 = 0; 
			double Le = sqrt(scni[i]->area);
			double c = sqrt(((lambda+2*mu)/rho));
			double P_b1 = b1*div_v*rho*Le*c;
			double eta = 0;
			double P_b2 = 0;
			if ( div_v < 0 ){
				eta = b1;
				P_b2 = Le*rho*(b2*b2*Le*div_v*div_v);
				eta -= b2*b2*Le*(1/c)*div_v;
			}





			double delta_t = (2.00/sqrt(((lambda + 2*mu)*MaxB)))*(sqrt(1+eta*eta)-eta);
			if ( delta_t <delta_t_min_i )
			{
				delta_t_min_i = delta_t;
			}



			if ((i == 124) && (call_count % 100 == 0)) {


			//m_foutput(stdout,stateNew[i]->W);
	
			stateNew[i]->F->me[2][1] = t_n_1;
			stateNew[i]->F->me[1][2] = stateNew[i]->critLambdaBar;
			stateNew[i]->F->me[2][0] = stateNew[i]->mu;
			stateNew[i]->F->me[0][2] = stateNew[i]->lambda;
			printf("lambda = %10.2E, mu = %10.2E \n",stateNew[i]->lambda,stateNew[i]->mu);
			//m_add(Sb_n_1,Savg_bond,Savg_bond);
				//m_add(Sc_n_1,Savg_conf,Savg_conf);
			snprintf(filename, 50, "strain_%d%s",print_count,".txt");
			mat2csv(stateNew[i]->F,"./History/Strain",filename);
			snprintf(filename, 50, "Bond_Stress_%d%s",print_count,".txt");
			mat2csv(stateNew[i]->Sb,"./History/Stress",filename);
			snprintf(filename, 50, "Conformational_Stress_%d%s",print_count,".txt");
			mat2csv(stateNew[i]->Sc,"./History/Stress",filename);
			int indx = 0;
			printf("max strain rate = %lf \n", v_max(stateNew[i]->lambdaDot,&indx));
			stateNew[i]->F->me[2][1] = 0;
			stateNew[i]->F->me[1][2] = 0;
			stateNew[i]->F->me[2][0] = 0;
			stateNew[i]->F->me[0][2] = 0;

			++print_count;


		}
			/* ------------------------------------------*/
			/* --------------Internal Force--------------*/
			/* ------------------------------------------*/


			/* Integration parameter */
			double intFactor = scni[i]->area;
			if( is_axi == 1){

				intFactor = intFactor * 2*PI*scni[i]->r;
			}

			/* Cauchy stress */
			// sigma->ve[0] = (Sc_n_1->me[0][0] + Sb_n_1->me[0][0] + mSigma + (P_b1+P_b2) )/1e6 ;
			// sigma->ve[1] = (Sc_n_1->me[1][1] + Sb_n_1->me[1][1] + mSigma + (P_b1+P_b2))/1e6 ;
			// sigma->ve[2] = (Sc_n_1->me[0][1] + Sb_n_1->me[0][1] )/1e6;
			// sigma->ve[3] = (Sb_n_1->me[2][2] + Sb_n_1->me[2][2] + mSigma + (P_b1+P_b2))/1e6 ;

			sigma->ve[0] = (stateNew[i]->sigma->me[0][0])/1e6 +(P_b1+P_b2);
			sigma->ve[1] = (stateNew[i]->sigma->me[1][1])/1e6 +(P_b1+P_b2);
			sigma->ve[2] = stateNew[i]->sigma->me[0][1]/1e6;
			sigma->ve[3] = (stateNew[i]->sigma->me[2][2])/1e6 +(P_b1+P_b2);

			// if ( i == 0)
			// {
			// 	m_foutput(stdout, stateNew[i]->Sb);
			// 	m_foutput(stdout, stateNew[i]->Sc);

			// }

			/*  Internal force vectors */
			gMat(G,stateNew[i]->invF,is_axi);
			mv_mlt(G,sigma,stressVoigt);
			sv_mlt(stateNew[i]->Jacobian,stressVoigt,stressVoigt);



		// push forward stress to new reference configuration
		if ( dim_piola == 4)
		{


			S11 = stressVoigt->ve[0] ; S12 = stressVoigt->ve[2];
			S21 = stressVoigt->ve[3]; S22 = stressVoigt->ve[1];

			F11 = scni[i]->F_r->me[0][0]; F12 = scni[i]->F_r->me[0][1];
			F21 = scni[i]->F_r->me[1][0]; F22 = scni[i]->F_r->me[1][1];

			stressVoigt->ve[0] = S11*F11 + S12*F12;
			stressVoigt->ve[1] = S21*F21 + S22*F22; 
			stressVoigt->ve[2] = S11*F21 + S12*F22;
			stressVoigt->ve[3] = S21*F11 + S22*F12;

		}else if ( dim_piola == 5)
		{


			S11 = stressVoigt->ve[0] ; S12 = stressVoigt->ve[2];
			S21 = stressVoigt->ve[3]; S22 = stressVoigt->ve[1];
			S33 = stressVoigt->ve[4] ;

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
	V_FREE(fIntTemp);
	M_FREE(G);
	V_FREE(sigma);
	V_FREE(stressVoigt);

	}  // end of parallel region 






	++call_count;

	return delta_t_min;
}
