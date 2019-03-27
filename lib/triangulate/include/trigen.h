
/*****************************************************************************/
/*                                                                           */
/*  (trigen.c)                                                              */
/*                                                                           */
/*  Program that calls triangle              */
/*                                                                           */
/*                                                                           */
/*  Takes in a boundary profile
	along with meshing options                             */
/*                                                                           */
/*****************************************************************************/
#ifndef TRIGEN_H_
#define TRIGEN_H_

#include "triangle.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define REAL double

typedef struct triangulateio * tri ;


// struct used to represent the triangle 
typedef struct _TRIANGLE{

	double * points;
	double * temperatures;
	int * boundary;
	int * triangles;
	int * pointmarkers;

	int num_points;
	int num_triangles;
	int num_boundary_points;



}TRIANGLE;

TRIANGLE  *  trigen(char * options, char * fileName);



#endif // end def TRIGEN_H_
