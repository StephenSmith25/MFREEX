
#include "../../MFREEX/include/Domain/Rigid/New_Rigid_Body.h"
#include <stdlib.h>
#include <math.h>

Rigid_Body * Create_Stretch_Rod(const double RADIUS, const double LENGTH, const int NUM_POINTS)
{

	const int DIM = 2;

	MAT * nodes = m_get(NUM_POINTS,DIM);

	MAT * srNodes_O = m_get(numPointsRod+1,2);
	for ( int i = 0 ; i < numPointsRod ; i++){
		double theta = -M_PI/2.00 + (M_PI/2/(numPointsRod-1))*i;
		srNodes->me[i][0] = stretchRodRad*cos(theta);
		srNodes->me[i][1] = 9.21+stretchRodRad*sin(theta);
		srNodes_O->me[i][0] = srNodes->me[i][0];
		srNodes_O->me[i][1] = srNodes->me[i][1];
	}
	srNodes->me[numPointsRod][0] = stretchRodRad;
	srNodes->me[numPointsRod][1] = 70;
	srNodes_O->me[numPointsRod][0] = stretchRodRad;
	srNodes_O->me[numPointsRod][1] = 70;


	Rigid_Body * body = malloc(1*sizeof(Rigid_Body));
	body->nodes = nodes;
	return body;
}


int Update_Stretch_Rod(Rigid_Body * sr, const double velocity,const double dt)
{


	return 0;

}