#ifndef GENERATE_MSCNI_H_
#define GENERATE_MSCNI_H_

#include "Integration/SCNI/generate_scni.h"
#include "matrix.h"
#include "matrix2.h"
#include <math.h>
#include "Integration/defgrad.h"
#define PI 3.14159265359


typedef struct MSCNI{

	double area; 
	IVEC * sfIndex ;
	VEC * phi;
	MAT * B ; 
	VEC * fInt;
	MAT * F_r;
	double * center;
	int index;

	double * stabalised_volumes;
	MAT ** stabalised_B;

}MSCNI;

typedef struct MSCNI_OBJ{

	MSCNI ** scni;
	int num_points;

}MSCNI_OBJ;


MSCNI_OBJ * generate_mscni(voronoi_diagram * voronoi, char * type, int is_stabalised, int is_AXI, meshfreeDomain * Mfree);

#endif