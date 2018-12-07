#ifndef GENERATE_VORONOI_H_
#define GENERATE_VORONOI_H_

#include "matrix.h"
#include "matrix2.h"


typedef struct voronoi_diagram                    
{
  MAT*                verticies;    

  int**               index;

  int* 				  num_cell_verticies;
  
  int 				  num_cells;


} voronoi_diagram;



voronoi_diagram * generate_voronoi(double *points, int * boundary, int numPoints, int numBoundary, int dim);
int print_voronoi_diagram(FILE *fp, voronoi_diagram * vor);

int free_voronoi_diagram(voronoi_diagram * vor);


#endif