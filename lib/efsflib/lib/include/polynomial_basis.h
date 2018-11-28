#ifndef POLYNOMIAL_BASIS_H_
#define POLYNOMIAL_BASIS_H_

#include "matrix.h"
#include "matrix2.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


int polynomial_basis(MAT * basis, double * point, int dim, char * order, int compute);

#endif