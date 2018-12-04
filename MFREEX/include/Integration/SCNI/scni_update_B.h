#ifndef SCNI_UPDATE_B
#define SCNI_UPDATE_B

#include "Integration/SCNI/generate_scni.h"
#include "matrix.h"
#include "matrix2.h"
#include "Integration/defgrad.h"
#define PI 3.14159265359

int scni_update_B(SCNI_OBJ * scni, VEC * disp, voronoi_diagram * voronoi, meshfreeDomain * Mfree, int is_AXI);

#endif