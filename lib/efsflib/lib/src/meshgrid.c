
#include "meshgrid.h"




MAT * meshgrid(double xmin,double xmax,int nx, double ymin, double ymax, int ny){

	MAT * xI = m_get(nx*ny,2);

	for ( int i = 0; i < nx ; ++i)
	{
		for ( int j = 0 ; j < ny ; ++j)
		{
			xI->me[i+i*(nx-1)+j][0] = xmin + (double)((xmax-xmin)/(nx-1))*i ;
			xI->me[i+i*(nx-1)+j][1] = ymin + (double)((ymax-ymin)/(ny-1))*j ;
		}
	}


	return xI;


}