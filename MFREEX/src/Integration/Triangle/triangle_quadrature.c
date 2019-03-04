
#include "Integration/Triangle/triangle_quadrature.h"



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




static inline int quadrature_points(double quad_array[][3], int quadorder){

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




// create quadrature points for each triangle 
QUAD_TRIANGLE * create_triangle_quadrature_points(double * points, int * triangles,
	int number_of_triangles, int * quad_orders )
{

	QUAD_TRIANGLE * quad_triangle = malloc(1*sizeof(QUAD_TRIANGLE));

	int n1, n2, n3;
	double x1,y1;
	double x2,y2;
	double x3,y3;
	double basis_functions[3][3];

	int num_quad_points_per_cell = 0;

	double cell_quad_points[20][3];


	int num_quad_points = 0;

	for ( int i = 0 ; i < number_of_triangles ; i++)
	{
		int quad_order = quad_orders[i];


		if ( quad_order == 1)
		{
			num_quad_points += 1;
		}else if ( quad_order == 2)
		{
			num_quad_points += 3;

		}else if (quad_order == 3)
		{
			num_quad_points += 4;

		}else if ( quad_order == 4)
		{
			num_quad_points += 6;

		}else{
			num_quad_points += 6;
		}
	}
	// estimate size
	MAT * quad_points = m_get(num_quad_points,2);
	VEC * quad_weights = v_get(num_quad_points);

	num_quad_points = 0; 



	for ( int i = 0 ; i < number_of_triangles ; i++)
	{
		// get 3 nodes of the triangle

		n1 = triangles[3*i]-1;
		n2 = triangles[3*i+1]-1;
		n3 = triangles[3*i+2]-1;


		int quad_order = quad_orders[i];

		int num_per_cell = quadrature_points(cell_quad_points, quad_order);
		// get coordinats of n1,n2 and n3

		x1 = points[2*n1];
		y1 = points[2*n1+1];
		x2 = points[2*n2];
		y2 = points[2*n2+1];
		x3 = points[2*n3];
		y3 = points[2*n3+1];




		for ( int k = 0 ; k < num_per_cell ; k++)
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





	}

	sv_mlt((1.00/2.000), quad_weights,quad_weights);
		
	quad_triangle->QUAD_POINTS = quad_points;
	quad_triangle->VOLUMES = quad_weights; 

	return quad_triangle;

}
