#ifndef TIMESTEP_H_
#define TIMESTEP_H_



#include "matrix.h"
#include "matrix2.h"
#include <stdlib.h>
#include <stdio.h>





typedef struct EXPLICIT_INTEGRATOR
{

	double dt;
	double t_n;
	double t_n_1;
	double t_max = 0;



}EXPLICIT_INTEGRATOR;

EXPLICIT_INTEGRATOR * NewExplicitIntegrator(void);

int InitialiseExplicitIntegration(EXPLICIT_INTEGRATOR * explicit_integrator);

int MakeTimeStep();



#endif