#ifndef PRESSURE_LOAD_H_
#define PRESSURE_LOAD_H_

#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "matrix2.h"
#include "mls_shapefunction.h"
#include "Boundary/amplitude.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

typedef enum PRESSURE_TYPE
{
	DIRECT_PRESSURE=1,
	CAVITY_PRESSURE=2,
}PRESSURE_TYPE;



/* Cavity pressure struct */
typedef struct CAVITYPRESSURE
{
	double A;
	double pLine; 
	double ChokedMassFlowRate;
	// lots more inputs 
	/*.
	.
	.
	.
	*/

}CAVITYPRESSURE;

/* Returns cavity pressure correponding to current parameters*/
double GetCavityPressure(CAVITYPRESSURE * cavity); // NOT Yet IMPLEMENTED
double GetChokedMassFlowRate(CAVITYPRESSURE * cavity); // NOT YET IMPLEMENTED 

typedef struct PRESSURE_LOAD
{

 	// Type of pressure load, e.g cavity or direct application
	PRESSURE_TYPE type;
	
	// else if cavity pressure, set upstream parameters
	CAVITYPRESSURE * cavity; 

	// if direct pressure find amplitude
	AMPLITUDE  amplitude;
	// if smooth
	SMOOTH_AMPLTIUDE * smooth_table;
	// if table
	TABLE * table;


	
	double pMag;


}PRESSURE_LOAD;




typedef struct pressure_boundary
{
	MAT * coords;
	IVEC * points;
	double * segment_weights;
	MAT * segment_normals;
	double magnitude;
	shape_function_container * sf_traction;
	int is_axi;

} pressure_boundary;

pressure_boundary * new_pressure_boundary(IVEC * points, meshfreeDomain * mFree);
int assemble_pressure_load(VEC * Fext, double Pressure, pressure_boundary * pB);
int update_pressure_boundary(pressure_boundary *pB, MAT * coords);
int free_pressure_boundary(pressure_boundary * pB);



#endif