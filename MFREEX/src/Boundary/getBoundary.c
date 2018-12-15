#include "Boundary/getBoundary.h"

static inline int findInt ( int * A, int val,int lArray){


	for ( int i = 0 ; i < lArray; i++){
		int temp = A[i];

		if ( temp == val){
			return i;
		}
	}




	return -1;
}

int getBoundary(IVEC ** bBoundary, int * Boundary,int numBoundary, int *pointmarkers, int numPoints, int nodeMarker )
{


	/*  Get the points which satsify the node marker */
	int tempArray[numPoints] ;
	int count = 0;
	for ( int i = 0 ; i < numPoints ; i++){
		int temp = pointmarkers[i];
		if ( temp == nodeMarker){
			tempArray[count] = i;
			count = count +1;
		}else{
			/*  Do nothing */
		}
	}
	/*  Initialise return vector */
	int numBnodes = count;
	IVEC * boundary = iv_get(numBnodes);

	/*  Find where the boundary starts, CW orientation */
	int bIndex = -1; 

	for ( int i = 0 ; i < numBnodes ; i++){
		int testNode = tempArray[i];

		int res = findInt(Boundary,testNode,numBoundary);
		int resm1 ;
		int resp1 ;
		if ( res == 0 ){
			resm1 = numBoundary;
			resp1 = res + 1;
		}
		else if( res == numBoundary -1){
			resp1 = 0;
			resm1 = res - 1;
		}else{
			resm1 = res -1;
			resp1 = res + 1;
		}

		if ( numBnodes < 3){

			if (( findInt(tempArray,Boundary[resm1],numBnodes) < 0) || 
					(findInt(tempArray,Boundary[resp1],numBnodes)> 0 )){
				bIndex = res;
				goto end;

			}
		}else{

			if (( findInt(tempArray,Boundary[resm1],numBnodes) < 0) && 
					(findInt(tempArray,Boundary[resp1],numBnodes)> 0 )){
				bIndex = res;
				goto end;
			}

		}
	}
end:;

	/*  Fill boundary vector */

	for ( int i = 0 ; i < numBnodes ; i++){
		int index = bIndex + i;
		if ( index > numBoundary - 1 ){
			index = index - numBoundary;	
		}
		boundary->ive[i] = Boundary[index] ;
	}

	*bBoundary = boundary;

	return 0;




	}
