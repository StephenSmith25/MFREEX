#ifndef MSCNI_UPDATE_B_H 
#define MSCNI_UPDATE_B_H

#include "Integration/SCNI/generate_mscni.h"
#include "matrix.h"
#include "matrix2.h"
#include "Integration/defgrad.h"
int mscni_update_B(MSCNI_OBJ * mscni_,  VEC * disp_, voronoi_diagram * vor, meshfreeDomain * Mfree, int is_AXI);


#endif