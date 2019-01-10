#include "Boundary/Displacement/enforceBC.h"






int enforceBC(EBC * ebc_, VEC * disp){


	// seperate disp into 1 and 2 directions
	VEC * disp_1 = v_get(disp->max_dim/2);
	VEC * disp_2 = v_get(disp->max_dim/2);


	for ( int i = 0 ; i < disp_1->max_dim ; i++){
		disp_1->ve[i] = disp->ve[2*i];
		disp_2->ve[i] = disp->ve[2*i+1];
	}

	VEC * uError = v_get(ebc_->nodes->max_dim);

	// 1 direction
	if ( ebc_->dofFixed == 1){

		mv_mltadd(ebc_->uBar1,disp_1,ebc_->phi,-1.00,uError);		
		mv_mlt(ebc_->P,uError,ebc_->uCorrect1);

		
		for ( int m = 0 ; m < ebc_->uCorrect1->max_dim ; m++){
			disp->ve[2*m] += ebc_->uCorrect1->ve[m];
		}

		// 2 direction 
	}else if(  ebc_->dofFixed == 2) {


		mv_mltadd(ebc_->uBar2,disp_2,ebc_->phi,-1.00,uError);	
		mv_mlt(ebc_->P,uError,ebc_->uCorrect2);
		for ( int m = 0 ; m < ebc_->uCorrect2->max_dim ; m++){
			disp->ve[2*m+1] += ebc_->uCorrect2->ve[m];
		}

	}else if ( ebc_->dofFixed == 3) {

		mv_mltadd(ebc_->uBar1,disp_1,ebc_->phi,-1.00,uError);		
		mv_mlt(ebc_->P,uError,ebc_->uCorrect1);

		for ( int m = 0 ; m < ebc_->uCorrect1->max_dim ; m++){
			disp->ve[2*m] += ebc_->uCorrect1->ve[m];
		}


		mv_mltadd(ebc_->uBar2,disp_2,ebc_->phi,-1.00,uError);	
		mv_mlt(ebc_->P,uError,ebc_->uCorrect2);
		for ( int m = 0 ; m < ebc_->uCorrect2->max_dim ; m++){
			disp->ve[2*m+1] += ebc_->uCorrect2->ve[m];
		}

		// Catch any errors
	}else{
		printf("DOF fixed variable not properly assisgned\n");
	}


	// Free any allocated memory that is no longer required

	v_free(disp_1);
	v_free(disp_2);
	v_free(uError);


	return 0;
}
