
#include "Boundary/iv_addNode.h"

int iv_addNode(IVEC * A, int a, char pos){






	int currentLength = A->max_dim;
	int newLength = currentLength + 1;
	IVEC * temp = iv_copy(A,IVNULL);
	iv_resize(A,newLength);
	if ( pos == 'e'){
		/*  Add at end */
		for ( int i = 0 ; i < newLength ; i++){
			if ( i < newLength - 1){
				A->ive[i] = temp->ive[i];
			}else if ( i == newLength - 1){
				A->ive[i] = a;
			}else{
				fprintf(stderr,"well something really went wrong here in iv_addNode \n");
			}
		}
	}else if (pos == 's'){
		/*  Add at start */	
		for ( int i = 0 ; i < newLength ; i++){
			if ( i== 0){
				A->ive[i] = a;
			}else if ( i > 0){
				A->ive[i] = temp->ive[i-1];
			}else{
				fprintf(stderr,"well something really went wrong here in iv_addNode \n");
			}
		}

	}else{
		/*  Char pos undefined */
		fprintf(stderr,"character position was not defined in iv_addNode\n");
	}

	IV_FREE(temp);



	return 0;
}
