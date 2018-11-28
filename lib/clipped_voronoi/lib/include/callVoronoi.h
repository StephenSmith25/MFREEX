#ifndef CALLVORONOI_H_
#define CALLVORONOI_H_



#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "gpc.h"

#define JCV_REAL_TYPE double
#define JCV_ATAN2 atan2
#define JCV_FLT_MAX 1.7976931348623157E+308
#include "jc_voronoi.h"

gpc_polygon ** callVoronoi(double points[], int num_points);






#endif
