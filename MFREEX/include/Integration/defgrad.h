#ifndef DEFGRAD_H_
#define DEFGRAD_H_

#include "matrix.h"
#include "matrix2.h"
#include "Integration/SCNI/generate_scni.h"

void get_defgrad(MAT * f,MAT * B,IVEC * neighbours, MAT * F_r, VEC * disp);


void get_dot_defgrad(MAT * f,MAT * B,IVEC * neighbours, MAT * F_r, VEC * velocity);




void get_grad_u(MAT * f, MAT * B, IVEC * neighbours, MAT * F_r, VEC * disp);







#endif