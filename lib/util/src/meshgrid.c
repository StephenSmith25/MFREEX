
#include "meshgrid.h"

#include "matrix.h"


MAT * meshgrid(double xmin,double xmax,int nx, double ymin, double ymax, int ny, int ** boundaryNodes, int * numBoundary){

	MAT * xI = m_get(nx*ny,2);

	int count = 0;

	*numBoundary = 2*nx + 2*(ny-2) ;
	int * boundaryNodes_out = malloc((2*nx + 2*(ny-2))*sizeof(int));



	for ( int i = 0; i < nx ; i++)
	{
		for ( int j = 0 ; j < ny ; ++j)
		{

			xI->me[i+i*(ny-1)+j][0] = (double)xmin + (double)((xmax-xmin)/(nx-1))*i ;
			xI->me[i+i*(ny-1)+j][1] = (double)ymin + (double)((ymax-ymin)/(ny-1))*j ;
		}
	}


	for ( int j = 0 ; j < ny ; j++)
	{

		int i = 0;
		int indx = i + i*(ny-1) + j;
		boundaryNodes_out[count] = i+i*(ny-1)+j;
		++count;

	}

	for ( int i = 1 ; i < nx ; i++)
	{

		int j = ny-1;
		int indx = i + i*(ny-1) + j;
		boundaryNodes_out[count] = i+i*(ny-1)+j;
		++count;


	}



	for ( int j = ny-2 ; j >0 ; j--)
	{

		int i = nx-1;
		int indx = i + i*(ny-1) + j;
		boundaryNodes_out[count] = i+i*(ny-1)+j;
		++count;

	}


	
	for ( int i = nx-1 ; i > 0 ; i--)
	{

		int j = 0;
		int indx = i + i*(ny-1) + j;

		boundaryNodes_out[count] = i+i*(ny-1)+j;
		++count;


	}

	*boundaryNodes = boundaryNodes_out;


	return xI;


}