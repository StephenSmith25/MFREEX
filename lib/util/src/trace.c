#include "trace.h"





double trace(MAT * matrix_){
	double returnVal = 0;

	if( matrix_->m!= matrix_->n){
		printf("finding trace of non square matrix, error\n");
		return 0;
	}else{
		returnVal = 0;
		for ( int i = 0 ; i < matrix_->m ; i++){
			
			returnVal += matrix_->me[i][i];
		}

		return returnVal;
	}
}
