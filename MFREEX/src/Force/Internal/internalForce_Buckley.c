#include "Force/Internal/internalForce_Buckley.h"


static int call_count;

double internalForce_ForceBuckley(VEC * Fint, SCNI_OBJ * scni_obj, VEC * disp, VEC * velocity,
	VEC * matParams,VEC * critLambdaParams, state_Buckley ** stateNew, state_Buckley ** stateOld,
	 int is_axi, int dim, double deltat)
{


	__zero__(Fint->ve,Fint->max_dim);



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

	
	// time step calculation
	double delta_t_min = 1000;



	// set number of threads
	omp_set_num_threads(1);

	int i;
#pragma omp parallel 
{


		// Internal force varaibles
		VEC * stressVoigt = v_get(dim_piola);
		VEC * sigma = v_get(dim_cauchy);
		VEC * fIntTemp = v_get(Fint->max_dim);

		// Gmat matrix
		MAT * G = m_get(dim_piola,dim_cauchy);

		// Time step
		double delta_t_min_i = 1000;

	
		MAT * B;
		IVEC * neighbours;
		MAT * F_r;
		MAT * Sc_n_1;
		MAT * Sb_n_1;


#pragma omp for nowait schedule(dynamic,4) 
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

			buckleyStress(stateNew[i],stateOld[i],matParams,critLambdaParams,deltat);

			/* ------------------------------------------*/
			/* --------------New time increment----------*/
			/* ------------------------------------------*/
			// m_add(Sb_n_1, Sc_n_1, cauchy_dev);


			// // Find hydrostatic and deviatoric stress

			// m_sub(cauchy_dev,stateOld[i]->dev_Stress,delta_cauchy_dev);

			// for ( int k = 0 ; k < dim_strain ; k++)
			// {
			// 	cauchy_hyd->me[i][i] = matParams->ve[5]*log(Jacobian_n_1);

			// }
			// m_sub(cauchy_hyd,stateOld[i]->hyd_Stress,delta_cauchy_hyd);
			//Find 1D frequency bounds
			// time step calculation
			// double B11 = 0;
			// double B22 = 0;
			// double MaxB = -1;
			// double rho = 1000e-9;


			// for ( int k = 0 ; k < (num_neighbours)*2 ; k++)
			// {	
			// 	B11 += (num_neighbours/rho)*(B->me[0][k]*B->me[0][k]);
			// 	B22 += (num_neighbours/rho)*(B->me[1][k]*B->me[1][k]);

			// }

			// MaxB = max(B11,B22);



			/* ------------------------------------------*/
			/* ---------------Bulk Damping---------------*/
			/* ------------------------------------------*/

			double div_v = stateNew[i]->div_v;
			double mSigma = stateNew[i]->mSigma;
			double Jacobian = stateNew[i]->Jacobian;


			double b1 = 0.06;
			double b2 = 1.44;
			double Le = 1.6e-3;
			double rho = 1380;
			//double c = sqrt(((lambda+2*mu)/rho));
			double c = 1400;
			double P_b1 = b1*div_v*rho*Le*c;
			double eta = b1;
			double P_b2 = 0;
			if ( div_v < 0 ){
				P_b2 = Le*rho*(b2*b2*Le*div_v*div_v);
				eta -= b2*b2*Le*(1/c)*div_v;
			}

			P_b2 = 0;
			P_b1 = 0;


			// double delta_t = (2.00/sqrt(((lambda + 2*mu)*MaxB)))*(sqrt(1+eta*eta)-eta);
			// if ( delta_t <delta_t_min_i )
			// {
			// 	delta_t_min_i = delta_t;
			// }


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

			sigma->ve[0] = stateNew[i]->sigma->me[0][0]/1e6;
			sigma->ve[1] = stateNew[i]->sigma->me[1][1]/1e6;
			sigma->ve[2] = stateNew[i]->sigma->me[0][1]/1e6;
			sigma->ve[3] = stateNew[i]->sigma->me[2][2]/1e6;



			/*  Internal force vectors */
			gMat(G,stateNew[i]->invF,is_axi);
			mv_mlt(G,sigma,stressVoigt);
			sv_mlt(stateNew[i]->Jacobian,stressVoigt,stressVoigt);

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
			// if ( delta_t_min_i < delta_t_min)
			// {
			// 	delta_t_min = delta_t_min_i;
			// }
		}

	// /*  Free allocated memory */
	// 	v_free(stressVoigt);
			V_FREE(fIntTemp);
			M_FREE(G);
			V_FREE(sigma);
			V_FREE(stressVoigt);

	// 	M_FREE(Fdot);

	}
	++call_count;

	return delta_t_min;
}
