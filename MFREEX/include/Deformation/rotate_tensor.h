#ifndef ROTATE_TENSOR_H_
#define ROTATE_TENSOR_H_

#include "matrix.h"
#include "matrix2.h"

int rotate_tensor(MAT * A, MAT * Q, MAT * temp,  MAT * OUT);
int un_rotate_tensor(MAT * A, MAT * Q, MAT * temp, MAT * OUT);


#endif