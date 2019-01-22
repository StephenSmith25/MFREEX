#ifndef MATERIAL_H_
#define MATERIAL_H_

#include "matrix.h"
#include "matrix2.h"


// stores the state of a material
typedef struct state
{

	MAT * F;
	MAT * D;


} state;


typedef struct state_Buckley
{
	// deformation gradient
	MAT * F;
	MAT * Fdot;
	MAT * Fn;

	MAT * invF;
	MAT * delta_F;
	MAT * delta_R;
	MAT * delta_U;
	MAT * delta_V;


	// isochoric deformation gradient
	MAT * Fbar;
	MAT * Fbardot;

	MAT * invFbar;

	// velocity gradinet
	MAT * L; 
	MAT * D;
	MAT * W;

	// Polar decomposition
	MAT * R;
	MAT * U;
	MAT * V;

	// Eigen values
	VEC * eigValDBar ;
	VEC * eigValVBar ;
	MAT * eigVecDBar  ;
	MAT * eigVecVBar ;

	// True strain
	MAT * true_strain;

	// Isochoric components
	MAT * Lbar; 
	MAT * Dbar;
	MAT * Wbar;

	// rotation
	MAT * Omega;

	// Deformation tensors
	MAT * Bbar;
	// bond and conformational stress
	MAT * Sb;
	MAT * Sc;
	MAT * sigma;

	// mean stress
	double mSigma;
	// stretch
	VEC * lambdaNBar;
	double lambdaNMax;
	// temperature
	double temperature;
	// critical network stretch
	double critLambdaBar;
	double Jacobian;
	double div_v; 

	MAT * h ; 
	MAT * GRAD_U ;

	// temp arrays
	MAT * temp;
	MAT * temp1;



} state_Buckley;

#endif