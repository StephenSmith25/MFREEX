#include "Force/Internal/internalForce_hyperelastic.h"


static int call_count;




int internalForce_hyperelastic(VEC * Fint, SCNI_OBJ * scni_obj, VEC * disp, VEC * matParams, char * Material, int is_axi, int dim){


	int i = 0 ;
	v_zero(Fint);

	// create pointer to correct material
	int (*mat_func_ptr)(VEC*,MAT*,VEC*) = NULL;

	if ( strcmp (Material, "SVK") == 0 )
	{
		mat_func_ptr = &SVK;

	}

	MAT * F = m_get(dim,dim);
	VEC * stressVoigt;

	// check if problem is axisymmetric
	if ( is_axi == 1){
		stressVoigt = v_get(5);
	}else{
		stressVoigt = v_get(dim*dim);
	}

	// loop over all integration points
	VEC * fIntTemp = v_get(Fint->max_dim);


	SCNI ** scni = scni_obj->scni;
	int num_int_points = scni_obj->num_points;


	for(i = 0 ; i < num_int_points ; i++){

		/*  Find deformation gradient */
		get_defgrad(F, scni[i]->B, disp, scni[i]->sfIndex);

		/* Integration parameter */
		double intFactor = scni[i]->area;
		if( is_axi == 1){

			intFactor = intFactor * 2*PI*scni[i]->center[0];
		}
		mat_func_ptr(stressVoigt,F,matParams);

		v_zero(scni[i]->fInt);
		vm_mlt(scni[i]->B,stressVoigt,scni[i]->fInt);

		for ( int k = 0 ; k < scni[i]->sfIndex->max_dim; k++){
			fIntTemp->ve[2*scni[i]->sfIndex->ive[k]] += intFactor * scni[i]->fInt->ve[2*k];
			fIntTemp->ve[2*scni[i]->sfIndex->ive[k]+1] += intFactor * scni[i]->fInt->ve[2*k+1];
		}
	}

	/*  Make this atomic or mutex so it is only done by one thread */
	v_add(fIntTemp,Fint,Fint);


	/*  Free allocated memory */
	v_free(stressVoigt);
	v_free(fIntTemp);
	M_FREE(F);

	++call_count;

	return 0;
}
