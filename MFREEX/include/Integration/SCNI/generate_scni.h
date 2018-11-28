#ifndef GENERATE_SCNI_H_
#define GENERATE_SCNI_H_

#include <stdlib.h>
#include "matrix.h"
#include "matrix2.h"
#include "mls_shapefunction.h"
#include "generate_voronoi.h"
#include "Bmat/generate_Bmat.h"
#include <math.h>
#include <omp.h>


typedef struct SCNI {
	double area; 
	IVEC * sfIndex ;
	VEC * phi;
	MAT * B ; 
	VEC * fInt;
	double * center;
	int index;
}SCNI;


SCNI ** generate_scni(voronoi_diagram * voronoi, char * type, int is_stabalised, int is_AXI, int dim, meshfreeDomain * Mfree);
int free_scni(SCNI * _scni);


#endif
