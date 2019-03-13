#ifndef CELLSEARCH_H_
#define CELLSEARCH_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"
#include "matrix2.h"

typedef struct NODELIST
{
	int node_number;

	struct NODELIST * next;

}NODELIST;

typedef struct CELL
{

	double * x;
	double * y;
	double * z;


	int num_nodes;


	NODELIST * nodes;


}CELL;


typedef struct CELLS
{
	CELL ** cells;
	int num_cells;
	int nx;
	int ny; 
}CELLS;

typedef struct BOUNDING_BOX
{
	double xmax,xmin;
	double ymax,ymin; 
	double zmax,zmin;



}BOUNDING_BOX;


/* Create bounding box to contain cells in */
BOUNDING_BOX * create_bounding_box(double xmin, double xmax, double ymin, double ymax,
	double zmin, double zmax);

/* Writes the bounding box to a file */
int  write_bounding_box_to_file();

/* Write cell to file */
int write_cells_to_file(CELLS * cells);

/*Create cells in a bounding box, each of a uniform size CELL_SIZE
and then place nodes in the cells */
CELLS * create_cells(BOUNDING_BOX * bounding_box, double * CELL_SIZE, int dim, MAT * NODES);


/* Iterate over each node in a cell*/
int iterate_over_cell(CELL * cell);


/* Add a particular node to a cell */
int add_node_to_cell(CELL * cell, int node_number);


/*Loop over each cell and adjust nodal cell location based on updated nodal coordinates*/
int update_cells(CELL ** cells, MAT * nodes);

#endif