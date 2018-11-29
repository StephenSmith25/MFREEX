#include "determinant.h"




double determinant ( MAT * matrix_){


	double returnVal = 0;  

	// check square matrix
	

	if ( matrix_->m != matrix_->n){
		printf("Trying to find determinant of non-square matrix \n");
		return 0;
	}
	else if ( matrix_->m == 2){


		returnVal = matrix_->me[0][0]*matrix_->me[1][1]  - matrix_->me[0][1]*matrix_->me[1][0];
		

		return returnVal;
	}

	else if ( matrix_->m == 3){


	returnVal =matrix_->me[0][0]* ((matrix_->me[1][1] * matrix_->me[2][2]) - (matrix_->me[1][2] * matrix_->me[2][1])) -
	matrix_->me[0][1]*((matrix_->me[1][0] * matrix_->me[2][2])-(matrix_->me[1][2] * matrix_->me[2][0])) +
	matrix_->me[0][2]* ((matrix_->me[1][0] * matrix_->me[2][1])- (matrix_->me[1][1] * matrix_->me[2][0]));

		return returnVal;
	}
	else{

		printf("trying to find determiannt of a matrix of an unsupported dimension \n");
		return 0; 
	}
	
	


}
