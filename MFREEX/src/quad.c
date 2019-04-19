#include "Integration/quad.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


static double QUAD_3D_P2_pts[][3] = {
	{0.5854101966249685,0.1381966011250105,0.1381966011250105},
	{0.1381966011250105,0.1381966011250105,0.1381966011250105},
	{0.1381966011250105,0.1381966011250105,0.5854101966249685},
	{0.1381966011250105,0.5854101966249685,0.1381966011250105}
};
static double QUAD_3D_P2_wts[] =
{
	0.250/6.00,0.250/6.00,0.250/6.00,0.250/6.00
};



static inline int lagrange_basis(double basis_functions[3][3], double xi, double eta)
{


    	//            (3)
		//            /\
		//           /  \
		//          /    \
		//         /      \
		//        /        \
		//       /          \
		//      /            \
		//     /              \
		//    /                \
		//   1----------------- 2




    //shape functions
    basis_functions[0][0]= 1.000 - xi - eta; //N1
    basis_functions[1][0]= xi; // N2
    basis_functions[2][0]=  eta; //N3

    // derivative with respect to xi
    basis_functions[0][1]= -1.000;
    basis_functions[1][1]= 1.000;
    basis_functions[2][1]= 0;

    // derivative with respect to eta
    basis_functions[0][2]= -1.000;
    basis_functions[1][2]=  0;
    basis_functions[2][2]= 1.000;

	return 0;
}




static inline int QUAD_2D_POINTS(double quad_array[][3], int quadorder){

	// TRIANGLE QUADRATURE
	// FROM DUNAVANT - HIGH DEGREE EFFICIENT SYMMETRICAL GAUSSIAN
	// QUADRATURE RULES FOR THE TRINAGLE, INT J. NUM. METH. ENG. 21(1985) 1129-1148

	int num_quad_points ;

	switch (quadorder) {

		case 1:
		num_quad_points = 1; 

		// x coords
		quad_array[0][0] = 0.33333333;

		// y coords
		quad_array[0][1] = 0.33333333;

		// weights
		quad_array[0][2] = 1.00000000;
		break;

		case 2:

		num_quad_points= 3; 
		// x coords
		quad_array[0][0] = 0.166666667;
		quad_array[1][0] = 0.166666667;
		quad_array[2][0] = 0.666666667;
				// y coords
		quad_array[0][1] = 0.166666667;
		quad_array[1][1] = 0.666666667;
		quad_array[2][1] = 0.166666667;

				// weights
		quad_array[0][2] = 0.33333333;
		quad_array[1][2] = 0.33333333;
		quad_array[2][2] = 0.33333333;

		break;

		case 3:

		num_quad_points= 4; 

		// x coords
		quad_array[0][0] = 0.333333333;
		quad_array[1][0] = 0.200000000;
		quad_array[2][0] = 0.200000000;
		quad_array[3][0] = 0.600000000;

		// y coords
		quad_array[0][1] = 0.333333333;
		quad_array[1][1] = 0.200000000;
		quad_array[2][1] = 0.600000000;
		quad_array[3][1] = 0.200000000;

		// weights
		quad_array[0][2] = -0.5625000;
		quad_array[1][2] = 0.52083333;
		quad_array[2][2] = 0.52083333;
		quad_array[3][2] = 0.52083333;

		break;

		case 4:


		num_quad_points= 6; 

		// x coords
		quad_array[0][0] = 0.44594849091597;
		quad_array[1][0] = 0.44594849091597;
		quad_array[2][0] = 0.10810301816807;
		quad_array[3][0] = 0.09157621350977;
		quad_array[4][0] = 0.09157621350977;
		quad_array[5][0] = 0.81684757298046;

		// y coords
		quad_array[0][1] = 0.44594849091597;
		quad_array[1][1] = 0.10810301816807;
		quad_array[2][1] = 0.44594849091597;
		quad_array[3][1] = 0.09157621350977;
		quad_array[4][1] = 0.81684757298046;
		quad_array[5][1] = 0.09157621350977;

				// weights
		quad_array[0][2] = 0.22338158967801;
		quad_array[1][2] = 0.22338158967801;
		quad_array[2][2] = 0.22338158967801;
		quad_array[3][2] = 0.10995174365532;
		quad_array[4][2] = 0.10995174365532;
		quad_array[5][2] = 0.10995174365532;

		break;

		default:


		printf(
			"Invalid quadrature order, maximum order = 4\ndefaulting to 4th order\n");
		

		num_quad_points= 6; 

		// x coords
		quad_array[0][0] = 0.44594849091597;
		quad_array[1][0] = 0.44594849091597;
		quad_array[2][0] = 0.10810301816807;
		quad_array[3][0] = 0.09157621350977;
		quad_array[4][0] = 0.09157621350977;
		quad_array[5][0] = 0.81684757298046;

				// y coords
		quad_array[0][1] = 0.44594849091597;
		quad_array[1][1] = 0.10810301816807;
		quad_array[2][1] = 0.44594849091597;
		quad_array[3][1] = 0.09157621350977;
		quad_array[4][1] = 0.81684757298046;
		quad_array[5][1] = 0.09157621350977;

				// weights
		quad_array[0][2] = 0.22338158967801;
		quad_array[1][2] = 0.22338158967801;
		quad_array[2][2] = 0.22338158967801;
		quad_array[3][2] = 0.10995174365532;
		quad_array[4][2] = 0.10995174365532;
		quad_array[5][2] = 0.10995174365532;

		break;

	}


	return num_quad_points;
}

// GETS THE QUADRATURE POINTS for a particular element 
QUAD * GetQuadPoints(MAT * nodes, ELEMENT * element, int order)
{
	QUAD * quad = NULL;

		if ( element->etype == TRIANGLES)
		{
			quad = QuadGetQuad1D(nodes, element->verticies, order);

		}else if ( element->etype == TETRAHEDRON)
		{
			quad = QuadGetQuad3D(nodes, element->verticies, order);
		}else if ( element->etype == LINE)
		{
			quad = QuadGetQuad1D(nodes, element->verticies, order);
		}else{
			fprintf(stderr,"Integration element not supported \n" );
		}


	return quad; 
}




QUAD * QuadGetQuad1D(MAT * nodes, int * verticies, int order)
{








	/* Returns 1D quad rule */

	return NULL; 

}





QUAD * QuadGetQuad2D(MAT * nodes, int * verticies, int order)
{


	// Initialise
	QUAD * quad = malloc(1*sizeof(QUAD));
	double quad_array[10][3];
	int num_quad_points = QUAD_2D_POINTS(quad_array,order);
	int n1, n2, n3;
	double x1,y1;
	double x2,y2;
	double x3,y3;
	double basis_functions[3][3];
	int num_quad_points_per_cell = 0;
	double cell_quad_points[20][3];

	// Set size of output matricies
	MAT * quad_points = m_get(num_quad_points,2);
	VEC * quad_weights = v_get(num_quad_points);

	// Vertex of triangls
	n1 = verticies[0];
	n2 = verticies[1];
	n3 = verticies[2];

	// get coordinats of n1,n2 and n3
	// X,Y coordinates of triangle vertiices
	x1 = nodes->me[n1][0];
	y1 = nodes->me[n1][1];
	x2 = nodes->me[n2][0];
	y2 = nodes->me[n2][1];
	x3 = nodes->me[n3][0];
	y3 = nodes->me[n3][1];

	// Place quadrature points in global reference frame 
	for ( int k = 0 ; k < num_quad_points ; k++)
	{
			lagrange_basis(basis_functions, cell_quad_points[k][0], cell_quad_points[k][1]);

			double j0[2][2];
			j0[0][0] = x1*basis_functions[0][1] + x2*basis_functions[1][1] + x3 * basis_functions[2][1];
			j0[0][1] = x1*basis_functions[0][2] + x2*basis_functions[1][2] + x3 * basis_functions[2][2];
			j0[1][0] = y1*basis_functions[0][1] + y2*basis_functions[1][1] + y3 * basis_functions[2][1];
			j0[1][1] = y1*basis_functions[0][2] + y2*basis_functions[1][2] + y3 * basis_functions[2][2];

			double detj0 = j0[0][0]*j0[1][1] - j0[0][1]*j0[1][0];


			quad_points->me[num_quad_points][0] = basis_functions[0][0]*x1 + basis_functions[1][0]*x2 
												+ basis_functions[2][0]*x3;
			quad_points->me[num_quad_points][1] = basis_functions[0][0]*y1 + basis_functions[1][0]*y2 
												+ basis_functions[2][0]*y3;		


			quad_weights->ve[num_quad_points] = cell_quad_points[k][2]*detj0;							 


			++num_quad_points;

	}
	sv_mlt((1.00/2.000), quad_weights,quad_weights);
		
	quad->points = quad_points;
	quad->weights = quad_weights; 
	quad->npoints = num_quad_points;

	return quad;

	return quad; 

}

//% Quadrature data for tetrahedron
//% Refs
//%  P Keast, Moderate degree tetrahedral quadrature formulas, CMAME 55: 339-348 (1986)
//%  O. C. Zienkiewicz, The Finite Element Method,  Sixth Edition,
QUAD * QuadGetQuad3D(MAT * nodes, int * verticies, int order)
{


	// NEED TO IMPLEMENT THIS 
	QUAD * quad = malloc(1*sizeof(QUAD));

	switch (order)
	{

		case (1): // N=4
		{
			quad->weights = malloc(4*sizeof(double));
			quad->points = malloc(3*4*sizeof(double));

			// place quadrature point the global frame

			// x 


			// y



			//


		}
		case (2):
		{

		}
		case (3):
		{



		}
		case (4):
		{

		}


	}
	/* Returns 1D quad rule */

	return NULL; 

}
