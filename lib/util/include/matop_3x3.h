#ifndef MATOP_3X3_H_
#define MATOP_3X3_H_

#include	"matrix.h"


MAT	*m_mlt_3x3(const MAT *A, const MAT *B, MAT *OUT);

MAT	*m_sub_3x3(const MAT *A, const MAT *B, MAT *OUT);
MAT	*m_add_3x3(const MAT *A, const MAT *B, MAT *OUT);

#endif