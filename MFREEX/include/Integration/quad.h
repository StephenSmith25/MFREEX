#ifndef QUAD_H_
#define QUAD_H_

#include "matrix.h"
#include "Mesh/Elements.h"


typedef struct QUAD_ {
	
	int order;
	int npoints;
	MAT * points;
	VEC * weights;

} QUAD;

void QuadFree(QUAD ** quad);

void QuadReset(void);


QUAD * QuadGetQuad1D(MAT * nodes, int * verticies, int order);
QUAD * QuadGetQuad2D(MAT * nodes, int * verticies, int order);
QUAD * QuadGetQuad3D(MAT * nodes, int * verticies, int order);
QUAD * GetQuadPoints(MAT * nodes, ELEMENT * element, int order);


#endif