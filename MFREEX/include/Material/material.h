#ifndef MATERIAL_H_
#define MATERIAL_H_

#include "matrix.h"
#include "matrix2.h"


// stores the state of a material in terms of its deformation, 
// and material properties



typedef struct state_variables
{


	/* ------------------------------------------*/
	/* -----------------Deformation--------------*/
	/* ------------------------------------------*/
	
	// Deformation gradient
	MAT * F;
	MAT * invF;

	// Determinant of F
	double Jacobian;

	// Deformation tensors
	MAT * B;
	MAT * C;

	// Polar Decomposition
	MAT * U;
	MAT * V;
	MAT * R;
	MAT * Vdot;



	// Rate measures
	MAT * L;
	MAT * D;
	MAT * W;
	MAT * Omega;
	double div_v;


	// Stress
	MAT * sigma;


	// rotated tensors
	MAT * d;
	MAT * T;
	MAT * omega; 



	// temperature
	double temperature;

	/* ------------------------------------------*/
	/* -----------------Plasticity--------------*/
	/* ------------------------------------------*/
	MAT * alpha;
	MAT * back_stress;
	MAT * eta;
	MAT * Deps;


	/* ------------------------------------------*/
	/* -----------------Buckley------------------*/
	/* ------------------------------------------*/

	// Eigen values
	VEC * eigValDBar;
	VEC * eigValVBar;
	MAT * eigVecDBar;
	MAT * eigVecVBar;
	VEC * lambdaDot;
	MAT * dbar;

	// Isochoric components
	MAT * Dbar;
	MAT * Lbar;
	MAT * Fbar;

	// Conformational branch network left cauchy green tensor
	MAT * Bbar;
	MAT * Dn;
	MAT * ep_n;
	// Network stretch
	VEC * lambdaNBar;


	// mean stress
	double mSigma;
	// max stretch
	double lambdaNMax;

	// currently conformational gamma
	double gamma; 

	// critical network stretch
	double critLambdaBar;
	// relaxation time
	double tau;

	// Buckley output stresses
	MAT * Sb;
	MAT * Sc;

	// rotated stresses
	MAT * Sb_R;
	MAT * Sc_R;

	/* ------------------------------------------*/
	/* -----------------Workspace----------------*/
	/* ------------------------------------------*/

	// Workspace matricies
	MAT * m_temp1;
	MAT * m_temp2;
	MAT * m_temp3;
	MAT * m_temp4;

	// Workspace vectors
	VEC * v_temp1;
	VEC * v_temp2;
	VEC * v_temp3;
	VEC * v_temp4;



} state_variables;

state_variables ** new_material_state(double * temperatures, int num_nodes, int is_buckley,
 int is_plastic, int dim, int is_AXI);

#endif