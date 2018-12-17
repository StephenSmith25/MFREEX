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
	MAT * invF;

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

	// True strain
	MAT * true_strain;

	// Isochoric components
	MAT * Fbar;
	MAT * invFBar;
	MAT * Lbar; 
	MAT * Dbar;
	MAT * Wbar;


	// Deformation tensors
	MAT * Bbar;
	// bond and conformational stress
	MAT * Sb;
	MAT * Sc;
	MAT * dev_Stress;
	MAT * hyd_Stress;
	// mean stress
	double mSigma;
	// stretch
	VEC * lambdaBar;
	double lambdaNMax;
	// temperature
	double temperature;
	// critical network stretch
	double critLambdaBar;
	double Jacobian;
	double div_v; 




} state_Buckley;

#endif