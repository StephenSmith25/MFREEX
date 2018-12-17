#include "Force/Internal/internalForce_Buckley.h"


static int call_count;

double internalForce_hyperelastic(VEC * Fint, SCNI_OBJ * scni_obj, VEC * disp, VEC * velocity,
	VEC * matParams, state_Buckley ** stateNew, state_Buckley ** stateOld, int is_axi, int dim, double deltat)
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


	double Kb = matParams->ve[9];

	// set number of threads
	omp_set_num_threads(8);


#pragma omp parallel 
{

		MAT * deltaF = m_get(dim_strain,dim_strain);
		MAT * invF_n_1 = m_get(dim_strain,dim_strain);
		MAT * invF_n_1_bar = m_get(dim_strain,dim_strain);
		MAT * Fdot = m_get(dim_strain,dim_strain);
		MAT * Fdot_bar = m_get(dim_strain,dim_strain);
		MAT * delta_R = m_get(dim_strain,dim_strain);
		MAT * delta_U = m_get(dim_strain,dim_strain);
		MAT * delta_V = m_get(dim_strain,dim_strain);
		MAT * delta_ep_true = m_get(dim_strain,dim_strain);
		MAT * temp = m_get(dim_strain,dim_strain);

		// Internal force varaibles
		VEC * stressVoigt = v_get(dim_piola);
		VEC * sigma = v_get(dim_cauchy);
		VEC * fIntTemp = v_get(Fint->max_dim);

		double delta_t_min_i = 1000;

		// Gmat matrix
		MAT * G = m_get(dim_piola,dim_cauchy);

		double Jacobian_n;
		double Jacobian_n_1;
		double mSigma;
		double div_v;

		MAT * B;
		IVEC * neighbours;
		MAT * F_r;



#pragma omp for nowait schedule(dynamic,4) 
		for(int i = 0 ; i < num_int_points ; i++){



			/* ------------------------------------------*/
			/* -----------------Deformation--------------*/
			/* ------------------------------------------*/

			// update material parameters
			// update materil parameters
			matParams->ve[11] = stateOld[i]->temperature+273.15;
			matParams->ve[12] = stateOld[i]->temperature+273.15;
			critLambdaParams->ve[5] = stateOld[i]->temperature + 273.15;		


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

			/* Find Fdot and Fbar dot */
			m_sub(stateNew[i]->F, stateOld[i]->F, delta_F);
			__smlt__(delta_F->base, 1.00/deltat, Fdot->base, dim_strain*dim_strain);

			// inverse deformation gradient
			m_inverse(F_n_1,invF_n_1);

			// polar decomposition to find delta R
			poldec(delta_F, delta_R, delta_U, delta_V);

			//------------------------//
			//      Velocity grad     //
			//------------------------//
			// Find velocity gradient 
			velocity_grad(stateNew[i]->L, stateNew[i]->D, stateNew[i]->W,Fdot,invF_n_1);
			// div_v
			if ( dim_strain == 2)
			{
				div_v = stateOld[i]->L->me[0][0] + stateOld[i]->L->me[1][1];

			}else if ( dim_strain == 3)
			{
				div_v = stateOld[i]->L->me[0][0] + stateOld[i]->L->me[1][1] + stateOld[i]->L->me[2][2];

			}


			// Get eigen values of polar decomposition
			tracecatch(
			symmeig(V,eigVecV,eigValV);,
			"internalForce");



			//------------------------//
			//   Strain increment     //
			//------------------------//


			// sm_mlt(deltat, stateOld[i]->D, delta_ep_true);
			// // update total strain
			// m_mlt(R, stateOld[i]->ep_true, temp);
			// mmtr_mlt(temp, R, stateOld[i]->ep_true);
			// m_add(stateOld[i]->ep_true, delta_ep_true, stateNew[i]->ep_true);


			// double delta_ep_vol = 0;

			// if ( dim == 2)
			// {
			// 	delta_ep_vol = (1.00/3.00)*(delta_ep_true->me[0][0] +delta_ep_true->me[1][1]);
			// 	delta_ep_true_vol->me[0][0] = delta_ep_vol;
			// 	delta_ep_true_vol->me[1][1] = delta_ep_vol;

			// }else if ( dim == 3)
			// {
			// 	delta_ep_vol = (1.00/3.00)*(delta_ep_true->me[0][0] +delta_ep_true->me[1][1] +delta_ep_true->me[2][2] );
			// 	delta_ep_true_vol->me[0][0] = delta_ep_vol;
			// 	delta_ep_true_vol->me[1][1] = delta_ep_vol;
			// 	delta_ep_true_vol->me[2][2] = delta_ep_vol;
			// }

			// m_sub(delta_ep_true,delta_ep_true_vol,delta_ep_true_iso);



			// Update Jacobian
			Jacobian_n = stateOld[i]->Jacobian;
			Jacobian_n_1 = Jacobian_n + Jacobian_n*div_v*deltat;

			/* ------------------------------------------*/
			/* ---------Isochoric   Deformation---------*/
			/* ------------------------------------------*/

			/* Distortional deformation */
			__smlt__(F_n_1->base, pow(Jacobian_n_1,-1.00/3.00), F_n_1_bar->base, dim_strain*dim_strain);
			// inv Fbar_n_1
			__smlt__(invF_n_1->base,pow(Jacobian_n_1,1.00/3.00), invF_n_1_bar->base, dim_strain*dim_strain);

			// F bar dot
			m_sub(F_n_1_bar, stateOld[i]->Fbar, Fdot_bar);
			__smlt__(Fdot->base, 1.00/deltat, Fdot_bar->base, dim_strain*dim_strain);

			velocity_grad(stateOld[i]->Lbar, stateOld[i]->Dbar, stateOld[i]->Wbar, Fdot_bar, invF_n_1_bar);

			// find eigen values of deviatoric deformation rate tensor
			symmeig(stateOld[i]->Dbar,eigVecD,eigValD);

			// update critical network stretch 
			stateOld[i]->critLambdaBar =lambdaCrit(stateOld[i]->critLambdaBar,eigValV, eigValD, critLambdaParams);


			/* ------------------------------------------*/
			/* ------------- --Find Stress--------------*/
			/* ------------------------------------------*/
			/*  Obtain stresses using explicit integration of stress rate */
			buckleyBond(stateNew[i],stateOld[i],matParams,deltat);
			// Conformational stress
			buckleyConf(stateNew[i],matParams,deltat);
			mSigma = log(Jacobian_n_1)*Kb;
			stateNew[i]->mSigma = mSigma;



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

			double b1 = 0.06;
			double b2 = 1.44;
			double Le = 1.6;
			double rho = 1380;
			//double c = sqrt(((lambda+2*mu)/rho));
			c = 1400;
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




			/* ------------------------------------------*/
			/* --------------Internal Force--------------*/
			/* ------------------------------------------*/


			/* Integration parameter */
			double intFactor = scni[i]->area;
			if( is_axi == 1){

				intFactor = intFactor * 2*PI*scni[i]->r;
			}


			/* Cauchy stress */
			sigma->ve[0] = (Sc_n_1->me[0][0] + Sb_n_1->me[0][0] + mSigma)/1e6 + qv;
			sigma->ve[1] = (Sc_n_1->me[1][1] + Sb_n_1->me[1][1] + mSigma)/1e6 + qv;
			sigma->ve[2] = (Sc_n_1->me[0][1] + Sb_n_1->me[0][1] )/1e6;
			sigma->ve[3] = (Sb_n_1->me[2][2] + Sb_n_1->me[2][2] + mSigma)/1e6 + qv;
			

			/*  Internal force vectors */
			gMat(G,invDefGrad);
			mv_mlt(G,sigma,stressVoigt);
			sv_mlt(Jacobian,stressVoigt,stressVoigt);

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
