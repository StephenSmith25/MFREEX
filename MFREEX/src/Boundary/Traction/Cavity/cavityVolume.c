#include "Boundary/Traction/Cavity/cavityVolume.h"


/*  Nodes is a integer vector of the nodal number that define the internal surface
 *  of the cavity. Coords is the current coords of these nodes*/

double cavityVolume(IVEC * nodes, MAT * coords){

	double volume = 0;

	for ( int i = 0 ; i < nodes->max_dim - 1 ; i++){
		// first point on segment
		double x1 = coords->me[nodes->ive[i]][0];
		double y1 = coords->me[nodes->ive[i]][1];
		//  second point on segment
		double x2 = coords->me[nodes->ive[i+1]][0];
		double y2 = coords->me[nodes->ive[i+1]][1];
		volume += fabs(PI*(1.000/3.000)*(x1*x1 + x1*x2 + x2*x2)*(y2-y1));
	}



	return volume ;
}
