#ifndef NEIGHBOURS_H_
#define NEIGHBOURS_H_

#include "matrix.h"
#include "matrix2.h"
#include <math.h>
#include "meshfree_structs.h"

IVEC * point_neighbours(double * x, meshfreeDomain * mfree);


IVEC * get_point_neighbours(IVEC * neighbours, double * x, meshfreeDomain * mfree);



#endif