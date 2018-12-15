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
	MAT * Fbar;
	double Jacobian;
	// velocity gradinet
	MAT * Dbar;
	MAT * Wbar;
	VEC * eigValDBar ;
	// Deformation tensors
	MAT * Bbar;
	// bond and conformational stress
	MAT * Sb;
	MAT * Sc;
	// mean stress
	double mSigma;
	// stretch
	VEC * lambdaBar;
	double lambdaNMax;
	// temperature
	double temperature;
	// critical network stretch
	double critLambdaBar;




} state_Buckley;

#endif