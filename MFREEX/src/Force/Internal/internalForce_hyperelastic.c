#include "Force/Internal/internalForce_hyperelastic.h"

static int call_count ;

int internalForce_hyperelastic(VEC * Fint, SCNI_OBJ * scni_obj, VEC * disp, VEC * matParams, char * Material, int is_axi, int dim){


	int i = 0 ;
	v_zero(Fint);


	__zero__(Fint->ve,Fint->max_dim);

	// create pointer to correct material
	int (*mat_func_ptr)(VEC*,MAT*,VEC*) = NULL;

	if ( strcmp (Material, "SVK") == 0 )
	{
		mat_func_ptr = &SVK;

	}

	int dim_v = 0;
	// check if problem is axisymmetric
	if ( is_axi == 1){
		dim_v = 5;
	}else{
		dim_v = dim*dim;
	}

	// loop over all integration points

	SCNI ** scni = scni_obj->scni;
	int num_int_points = scni_obj->num_points;


	omp_set_num_threads(3);

#pragma omp parallel 
{

	MAT * F = m_get(dim,dim);
	VEC * stressVoigt = v_get(dim_v);
	VEC * fIntTemp = v_get(Fint->max_dim);

#pragma omp for nowait schedule(dynamic,2)
	for(i = 0 ; i < num_int_points ; i++){

		/*  Find deformation gradient */
		get_defgrad(F, scni[i], disp);

		/* Integration parameter */
		double intFactor = scni[i]->area;
		if( is_axi == 1){

			intFactor = intFactor * 2*PI*scni[i]->center[0];
		}
		if ( call_count >49999 )
		{
			if ( i == 0){
			printf("F_r ");
			m_foutput(stdout, scni[i]->F_r);

			printf("F ");
			m_foutput(stdout, F);
			}
		}

		mat_func_ptr(stressVoigt,F,matParams);
		__zero__(scni[i]->fInt->ve,scni[i]->fInt->max_dim);

		vm_mlt(scni[i]->B,stressVoigt,scni[i]->fInt);

		for ( int k = 0 ; k < scni[i]->sfIndex->max_dim; k++){
			fIntTemp->ve[2*scni[i]->sfIndex->ive[k]] += intFactor * scni[i]->fInt->ve[2*k];
			fIntTemp->ve[2*scni[i]->sfIndex->ive[k]+1] += intFactor * scni[i]->fInt->ve[2*k+1];
		}
	}

	/*  Make this atomic or mutex so it is only done by one thread */
	#pragma omp critical
		__add__(fIntTemp->ve, Fint->ve, Fint->ve, Fint->max_dim);

	/*  Free allocated memory */
	v_free(stressVoigt);
	v_free(fIntTemp);
	M_FREE(F);

}
	++call_count;

	return 0;
}
