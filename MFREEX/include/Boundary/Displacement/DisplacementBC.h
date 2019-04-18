#ifndef DISPLACEMENTBC_H_
#define DISPLACEMENTBC_H_

#include "Boundary/amplitude.h"

typedef enum BC_TYPE
{
	DOF_FIXED=1,
	DOF_PRESCRIBED=2,
}BC_TYPE;

typedef enum DOF_TYPE
{
	DISPLACEMENT=1,
	VELOCITY=2,
}DOF_TYPE;


typedef enum DIRECTION
{
	X=1,
	Y=2,
	Z=3
}DIRECTION;

typedef struct DOF_CONSTRAINT
{
	DIRECTION dir;
	BC_TYPE type;

	// how the prescribed value is 
	AMPLITUDE amplitude;

	// If instant it will be prescribed from t=0;
	double value;
	SMOOTH_AMPLTIUDE * smooth_amplitude;
	TABLE * table_amplitude;

	// pointer to next dof constrained 
	struct DOF_CONSTRAINT * next;



}DOF_CONSTRAINT;

void PrintConstraintType(DOF_CONSTRAINT * dof_constraint);

#endif