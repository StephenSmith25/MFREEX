
#include "Boundary/Contact/contactDetection.h"


double contactDetection(MAT * point, MAT * masterNodes,MAT  * msNormal){




	double distMin = 1e4;
	int closestNode = -1;
	double segNormal[2];
	double segLength = 0;
	int isContact = -1;
	double delta ;

	/*  Have to test each for the closest segment */
	/*  Find closest master nodes*/
	for ( int i = 0 ; i < masterNodes->m ; i++){
		/*  Find closest master nodes */

		double dist = sqrt(pow(masterNodes->me[i][0]-point->me[0][0],2) + pow(masterNodes->me[i][1]-point->me[0][1],2));


		if ( dist < distMin){

			/*  Find minimum distance */
			closestNode = i;
			distMin = dist;
		}



	}
	/*  Now find out which segment the closest point projection lies in */


	int closestSeg = -1;

	/*  Test the dot product to the left and right of the closest node */

	if ( closestNode == 0){
		closestSeg = closestNode + 1;
		segNormal[0] = 0;
		segNormal[1] = -1;

		double eA[2] = {masterNodes->me[closestNode+1][0]-masterNodes->me[closestNode][0],
			masterNodes->me[closestNode+1][1]-masterNodes->me[closestNode][1]};
		double masterDist[2] = {point->me[0][0] - masterNodes->me[closestNode][0],
			point->me[0][1] - masterNodes->me[closestNode][1]};
		double dot2 = masterDist[0]*eA[0] + masterDist[1]*eA[1];


		if ( dot2 < 0){
			closestSeg = closestNode;
		}

	}else if ( closestNode == (masterNodes->m -1)){
		closestSeg = closestNode - 1;
		segNormal[0] = 1;
		segNormal[1] = 0;	

		double eAm1[2] = {masterNodes->me[closestNode][0]-masterNodes->me[closestNode-1][0],
			masterNodes->me[closestNode][1]-masterNodes->me[closestNode-1][1]};
		double masterDist[2] = {point->me[0][0] - masterNodes->me[closestNode][0],
			point->me[0][1] - masterNodes->me[closestNode][1]};

		double dot1 = masterDist[0]*eAm1[0] + masterDist[1]*eAm1[1];

		if ( dot1 > 0){
			closestSeg = closestNode;
		}

	}else{
		/* Find tangent vector to segment A-1---A  */
		double eAm1[2] = {masterNodes->me[closestNode][0]-masterNodes->me[closestNode-1][0],
			masterNodes->me[closestNode][1]-masterNodes->me[closestNode-1][1]};
		/* Find tangent vector to segment A---A+1  */
		double eA[2] = {masterNodes->me[closestNode+1][0]-masterNodes->me[closestNode][0],
			masterNodes->me[closestNode+1][1]-masterNodes->me[closestNode][1]};
		/*  Can check both sides */

		double masterDist[2] = {point->me[0][0] - masterNodes->me[closestNode][0],
			point->me[0][1] - masterNodes->me[closestNode][1]};


		double dot1 = masterDist[0]*eAm1[0] + masterDist[1]*eAm1[1];
		double dot2 = masterDist[0]*eA[0] + masterDist[1]*eA[1];

		if (( dot1 <= 0) && (dot2 < 0)){
			closestSeg = closestNode - 1;	

			segLength = sqrt(pow(eAm1[0],2) + pow(eAm1[1],2));
			segNormal[0] = eAm1[1]/segLength;
			segNormal[1] = -eAm1[0]/segLength;

		}
		if (( dot1 > 0) && (dot2 >= 0)){
			closestSeg = closestNode +1;

			segLength = sqrt(pow(eA[0],2) + pow(eA[1],2));
			segNormal[0] = eA[1]/segLength;
			segNormal[1] = -eA[0]/segLength;
		}
		if (( dot1 <= 0) && (dot2 >= 0)){
			closestSeg = closestNode +1;
			segLength = sqrt(pow(eA[0],2) + pow(eA[1],2));
			segNormal[0] = eA[1]/segLength;
			segNormal[1] = -eA[0]/segLength;

		}
		if (( dot1 > 0) && (dot2 < 0)){
			closestSeg = closestNode ;

			segLength = sqrt(pow(eA[0],2) + pow(eA[1],2));
			segNormal[0] = eA[1]/segLength;
			segNormal[1] = -eA[0]/segLength;
		}


	}

	/*  Now closest distance */
	/*  Find xi  */
	/*  Segment is represented by 
	 *
	 *	y = y0 + xi ( a )
	 *  a = y1 - y0
	 *
	 *  closest point is defined by
	 *
	 *	(x-y). (y),xi  = 0 defines a normal 
	 *
	 *  (x - y0 + xi(a)) . a = 0
	 *
	 *  x - y0 = b;
	 *
	 *
	 *  (b + xi(a)) . a = 0
	 *
	 *  b.a = -xi ( a. a )
	 *
	 *  */

	double closestDistance = 0;

	if ( closestSeg != closestNode){
		/*  Parametric representation of segment */
		double a[2] = {masterNodes->me[closestSeg][0]-masterNodes->me[closestNode][0],
			masterNodes->me[closestSeg][1]-masterNodes->me[closestNode][1]};
		double b[2] = {point->me[0][0] - masterNodes->me[closestNode][0],
			point->me[0][1] - masterNodes->me[closestNode][1]};
		double aDotA = a[0]*a[0] + a[1]*a[1];
		double bDotA = b[0]*a[0] + b[1]*a[1];
		double xi = ( bDotA / aDotA);

		double closestPoint[2] = {masterNodes->me[closestNode][0] + xi*a[0],
			masterNodes->me[closestNode][1] + xi*a[1]};
		closestDistance = sqrt(   pow(closestPoint[0] - point->me[0][0],2)   + pow(closestPoint[1] - point->me[0][1],2)     );


		double vec2slave[2] = {closestPoint[0]-point->me[0][0],closestPoint[1]-point->me[0][1]};
		delta = vec2slave[0]*segNormal[0] + vec2slave[1]*segNormal[1];





	}else{
		double closestPoint[2] = {masterNodes->me[closestNode][0],masterNodes->me[closestNode][1]};
		closestDistance = sqrt(   pow(closestPoint[0] - point->me[0][0],2)   + pow(closestPoint[1] - point->me[0][1],2)     );

		double vec2slave[2] = {closestPoint[0]-point->me[0][0],closestPoint[1]-point->me[0][1]};
		delta = vec2slave[0]*segNormal[0] + vec2slave[1]*segNormal[1];
	}


	/*  point lies between closestSegNode and closestNode */

	//printf("segNormal = {%lf,%lf}\n",segNormal[0],segNormal[1]);
	/*  Return closest segment, and if it was in contact? */
	//	delta = max(delta,0);

	msNormal->me[0][0] = segNormal[0];
	msNormal->me[0][1] = segNormal[1];
	return delta;



	/*  Print return the gap g.n  */
}
