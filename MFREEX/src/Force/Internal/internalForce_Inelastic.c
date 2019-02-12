#include "Force/Internal/internalForce_Inelastic.h"
#include "Deformation/update_Polar_Decomposition.h"
#include "Deformation/rotate_tensor.h"
#include "Material/Buckley/buckleyStress.h"


static int call_count;
static int print_count = 1;

double 
internalForce_Inelastic(VEC * Fint, SCNI_OBJ * scni_obj,
	VEC * disp, VEC * velocity,
	VEC * matParams, 
	state_variables ** stateNew, 
	state_variables ** stateOld,
	int is_axi, int dim, 
	double DT, double t_n_1, 
	char * Material)
{



		// zero internal force input
	__zero__(Fint->ve,Fint->max_dim);
	double rho;

	// create pointer to correct material function
	int (*mat_func_ptr)(state_variables *,state_variables * ,VEC *, double) = NULL;

	if ( strcmp (Material, "BUCKLEY") == 0 )
	{
		mat_func_ptr = &buckleyStress;
		rho = matParams->ve[31];

	}else if( strcmp(Material, "J2") == 0 )
	{
		//mat_func_ptr = &yeoh;


	}

		// Get dimensions of working vectors and matricies
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
	omp_set_num_threads(16);

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



			B = scni[i]->B;
			neighbours = scni[i]->sfIndex;
			F_r = scni[i]->F_r;
			int num_neighbours = neighbours->max_dim;

			/* ------------------------------------------*/
			/* -----------------Deformation--------------*/
			/* ------------------------------------------*/
			
			/*  Find deformation gradient */
			get_defgrad(stateNew[i]->F, B, neighbours,F_r,disp);

			/*  Find Rate of deformation tensors at n+h*/
			velocity_grad(stateNew[i],stateOld[i],DT,0.500);

			// Find inverse deformation gradient at n+1
			m_inverse_small(stateNew[i]->F, stateNew[i]->invF);

			// Find Jacobian at n+1
			stateNew[i]->Jacobian = determinant(stateNew[i]->F);


			//------------------------------------------//
			//        Update Polar Decomposition        //
			//------------------------------------------//
			update_Polar_Decomposition(stateNew[i], stateOld[i], DT);

			// remove rotation from D to get d
			un_rotate_tensor(stateNew[i]->D, stateNew[i]->R, 
				stateNew[i]->m_temp1, stateNew[i]->d);


			//--------------------------------------------//
			//         Stress from constitutive law       //
			//--------------------------------------------//

			mat_func_ptr(stateNew[i],stateOld[i],matParams,DT);



			/* ------------------------------------------*/
			/* --------------Print Outputs--------------*/
			/* ------------------------------------------*/


			if ((i == 182) && (call_count % 100 == 0)) {


				stateNew[i]->F->me[2][1] = t_n_1;

				snprintf(filename, 50, "strain_%d%s",print_count,".txt");
				mat2csv(stateNew[i]->F,"./History/Strain",filename);
				//snprintf(filename, 50, "Stress_%d%s",print_count,".txt");
				//mat2csv(stateNew[i]->sigma,"./History/Stress",filename);

				snprintf(filename, 50, "Bond_Stress_%d%s",print_count,".txt");
				mat2csv(stateNew[i]->Sb,"./History/Stress",filename);
				snprintf(filename, 50, "Conformational_Stress_%d%s",print_count,".txt");
				mat2csv(stateNew[i]->Sc,"./History/Stress",filename);
				stateNew[i]->F->me[2][1] = 0;
				++print_count;


			}


			/* ------------------------------------------*/
			/* -----------------Damping -----------------*/
			/* ------------------------------------------*/

			double b1 = 0.01;
			double b2 = 1.2;
			double Le = 2.42;
			double Cd = 1400;
			double div_v = stateNew[i]->div_v;
			double qv =  rho*Le*b1*Cd * div_v;
			if ( div_v < 0)
			{
				qv += rho*Le*(b2 * (Le/1000) * pow(div_v,2)) ;
			}





			/* ------------------------------------------*/
			/* --------------Internal Force--------------*/
			/* ------------------------------------------*/


			/* Integration parameter */
			double intFactor = scni[i]->area;
			if( is_axi == 1){

				intFactor = intFactor * 2*PI*scni[i]->r;
			}


			//-------------------------------------------------//
			//         Get Voigt Piola Kirchoff stress         //
			//-------------------------------------------------//


			sigma->ve[0] = (stateNew[i]->sigma->me[0][0])/1e6 + qv;
			sigma->ve[1] = (stateNew[i]->sigma->me[1][1])/1e6 + qv;
			sigma->ve[2] = stateNew[i]->sigma->me[0][1]/1e6;
			sigma->ve[3] = (stateNew[i]->sigma->me[2][2])/1e6 + qv;



			/*  Internal force vectors */
			gMat(G,stateNew[i]->invF,is_axi);
			mv_mlt(G,sigma,stressVoigt);
			sv_mlt(stateNew[i]->Jacobian,stressVoigt,stressVoigt);

		
			//-------------------------------------------------//
			//            Assemble Internal Force              //
			//-------------------------------------------------//

			vm_mlt(scni[i]->B,stressVoigt,scni[i]->fInt);
			for ( int k = 0 ; k < num_neighbours; k++){
				fIntTemp->ve[2*scni[i]->sfIndex->ive[k]] += intFactor * scni[i]->fInt->ve[2*k];
				fIntTemp->ve[2*scni[i]->sfIndex->ive[k]+1] += intFactor * scni[i]->fInt->ve[2*k+1];
			}


			//-------------------------------------------------//
			//                  State updates                  //
			//-------------------------------------------------//

			// n+1 --->>>> n 
			m_copy(stateNew[i]->F,stateOld[i]->F);
			m_copy(stateNew[i]->V,stateOld[i]->V);
			m_copy(stateNew[i]->R,stateOld[i]->R);
			m_copy(stateNew[i]->invF,stateOld[i]->invF);			
			stateOld[i]->Jacobian = stateNew[i]->Jacobian;
			stateOld[i]->div_v = stateNew[i]->div_v;

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



	return 0.000;
}

