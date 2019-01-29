#include "internalForce.h"

#include "SVK.h"






int internalForce(VEC * Fint, SCNI * scni, VEC * disp, VEC * matParams, int numnode ){


	int i = 0 ;

	v_zero(Fint);
	VEC * stressVoigt;
	/*  If return type = C; have to find G and J, and multiply cauchy stress voigt by it */
	if ( AXI == 1){
		stressVoigt = v_get(5);
	}else{
		stressVoigt = v_get(4);
	}
	VEC * fIntTemp = v_get(2*numnode);
	for(i = 0 ; i < numnode ; i++){
		/*  Find deformation gradient */
		smoothedDefGrad(&scni[i], disp, numnode);
		double intFactor = scni[i].area;	
		if( AXI == 1){

			intFactor = intFactor * 2*PI*scni[i].center[0];
		}
		mooneyRivlin(stressVoigt,scni[i].F,matParams);
		v_zero(scni[i].fInt);
		vm_mlt(scni[i].B,stressVoigt,scni[i].fInt);

		for ( int k = 0 ; k < scni[i].sfIndex->max_dim; k++){
			fIntTemp->ve[2*scni[i].sfIndex->ive[k]] += intFactor * scni[i].fInt->ve[2*k];
			fIntTemp->ve[2*scni[i].sfIndex->ive[k]+1] += intFactor * scni[i].fInt->ve[2*k+1];
		}
	}

	/*  Make this atomic or mutex so it is only done by one thread */
	v_add(fIntTemp,Fint,Fint);

	/*  Free allocated memory */
	v_free(stressVoigt);
	v_free(fIntTemp);

	return 0;
}
