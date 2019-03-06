#ifndef MESHFREE_STRUCTS_H_
#define MESHFREE_STRUCTS_H_

#include "matrix.h"
#include "matrix2.h"


enum SUPPORT_TYPE {RADIAL,RECTANGULAR,ELLIPTICAL};



typedef struct meshfreeDomain
{
	MAT * nodes;
	VEC * di;
	MAT * di_tensor;

	MAT ** MI;

	double dmax_radial ;
	double * dmax_tensor;
	double beta ;

	char * kernel_shape;
	char * basis_type;
	char * weight_function;


	enum SUPPORT_TYPE kernel_support;
	int num_nodes;
	int dim;
	int IS_AXI;
	int is_constant_support_size;
} meshfreeDomain;

typedef struct shape_function{


	// construction point
	double * x;
	IVEC * neighbours;
	VEC * phi;
	MAT * dphi;
	MAT * d2phi;
	int compute;
	int dim;  


	/* INTERNAL VARIABLES */ 

	// initialise shape function calculation matricies (these will be free'd
	MAT * A ;
	MAT * B ;
	MAT * basis_xi ;
	VEC * p_xi ;
	MAT * ppT ;
	VEC * bi ;
	VEC * weights ;
	MAT * invA ;
	VEC * inter ;
	MAT * basis ;
	VEC * p ;
	PERM * pivot ;
	MAT * LU_A ;

	// DERIVATIVES
	VEC * dphi_dk;
	VEC * dp_dk ;
	MAT ** dA_dk ;
	MAT ** dB_dk ;

	// get second derivatives of phi
	MAT ** d2A_dk_dj ;
	MAT ** d2B_dk_dj ;
	VEC * d2phi_dk_dj;

	// TEMPORY
	VEC * v_inter_1 ;
	VEC * v_inter_2 ;
	VEC * v_inter_3 ;



	// LU decompostion variables
	VEC * gamma ;
	VEC * gamma_k ;
	VEC * gamma_kj ;
	MAT * gamma_k_m ;
	VEC * RHS ;


	/* END INTERNAL VARIABLES */



} shape_function;

typedef struct shape_function_container{

	int number_of_points;
	int compute;
	shape_function ** sf_list;

} shape_function_container;
#endif