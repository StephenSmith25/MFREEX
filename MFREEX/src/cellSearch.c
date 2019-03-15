#include "cellSearch.h"
#include <stdbool.h>


static inline void deleteList(NODELIST * prev)
{
	// save the node we wish to delete in a temporary pointer
	NODELIST * temp = prev->next;

	// set the previous nodes next pointer to point to the node after the node we wish to delete
	prev->next = temp->next;

}
static inline void insertList(CELL * cell, NODELIST * newpoint)
{

	// link the new item to point ot the head of the list
	newpoint->next = cell->nodes;

	// set of head of the list be our new item
	cell->nodes = newpoint;


}


CELLS * create_cells(BOUNDING_BOX * bounding_box, double * CELL_SIZE, int dim, MAT * NODES)
{



	CELLS * cells = malloc(1*sizeof(CELLS));


	if ( dim == 2)
	{


		// GET CELL DIMENSIONS
		double CX = CELL_SIZE[0];
		double CY = CELL_SIZE[1];




		// Find nx
		double Bx_min = bounding_box->xmin;
		double Bx_max = bounding_box->xmax;

		int nx = ceil(fabs(Bx_max - Bx_min) / CX );

		// Find ny
		double By_min = bounding_box->ymin;
		double By_max = bounding_box->ymax;


		int ny = ceil(fabs(By_max - By_min) / CY );

		cells->nx = nx;
		cells->ny = ny;

		// Generate the nx and ny cells
		cells->num_cells = ny*nx;
		cells->cells = malloc(nx*sizeof(CELL*));
		for ( int i =0 ; i < nx ; i++)
		{
			cells->cells[i] = malloc(ny*sizeof(CELL));
		}


		for ( int i = 0 ; i < nx ; i++ )
		{
			for (int j = 0 ; j < ny ; j++)
			{
				cells->cells[i][j].x = malloc(2*sizeof(double));
				cells->cells[i][j].y = malloc(2*sizeof(double));


				/*Create cell limits*/
				cells->cells[i][j].x[0] = i*CX;
				cells->cells[i][j].x[1] = (i+1)*CX;

				cells->cells[i][j].y[0] = j*CY;
				cells->cells[i][j].y[1] = (j+1)*CY;

				// CREATE CELL NODE LIST
				cells->cells[i][j].nodes = NULL;

				cells->cells[i][j].num_nodes = 0;


				for ( int k = 0 ; k < NODES->m ; k++)
				{
					double x_test = NODES->me[k][0];
					double y_test = NODES->me[k][1];

					if (( x_test >= cells->cells[i][j].x[0]) && ( x_test <= cells->cells[i][j].x[1]))
					{
						if (( y_test >= cells->cells[i][j].y[0]) && ( y_test <= cells->cells[i][j].y[1]))
						{
							/*Create nodelist object */
							NODELIST * node_i = (NODELIST * )malloc(1*sizeof(NODELIST));
							node_i->node_number = k;

							// NODE K IS INSIDE THIS CELL
							if ( cells->cells[i][j].num_nodes == 0)
							{
								cells->cells[i][j].nodes = node_i;
								cells->cells[i][j].nodes->next = NULL;
							}else{
								node_i->next = cells->cells[i][j].nodes;
								cells->cells[i][j].nodes = node_i;
							}
							
							++cells->cells[i][j].num_nodes ;


							/*Remove that node from possible search candidates  
							- to avoid duplicated nodes*/
							NODES->me[k][0] = -100000000;
							NODES->me[k][1] = -100000000;


						}

					}

				}


			}
		}




	}else if ( dim ==3 )
	{

		// Find nx
		// Find ny
		// Find nz

	}



	 return cells;
}


/* Create bounding box to contain cells in */
BOUNDING_BOX * create_bounding_box(double xmin, double xmax, double ymin, double ymax,
	double zmin, double zmax)
{



	BOUNDING_BOX * bounding_box = malloc(1*sizeof(BOUNDING_BOX));


	bounding_box->xmin = xmin; bounding_box->xmax = xmax;
	bounding_box->ymin = ymin; bounding_box->ymax = ymax;
	bounding_box->zmin = zmin; bounding_box->zmax = zmax;



	return bounding_box;
}

static inline void deleteCell()
{


}

static inline void addCell(int cell_number)
{


}

active_cell * get_active_cells(CELLS * cells, int * num_active_cells)
{

	int i,j,k; 

	int count_num_active_cells = 0;
	int nx = cells->nx;
	int ny = cells->ny;

	active_cell * active_cells = NULL;
	active_cell * previous = NULL;

	for ( int i = 0 ; i < cells->nx ; i++)
	{

		for ( int j = 0 ; j < cells->ny ; j++)
		{
			if ( cells->cells[i][j].nodes != NULL)
			{
				// add this cell to the active cells list
				active_cell * active_cell_i = malloc(1*sizeof(active_cell));
				active_cell_i->cell_number = (i)*nx + j;

				if (previous==NULL)
				{
					active_cells = active_cell_i;
				}else{
					previous->next = active_cell_i;
					previous->next->next = NULL;
				}	

				previous = active_cell_i;

				++count_num_active_cells;
			}
		}

	}


	// returns a sorted array of active cell indicies
	*num_active_cells = count_num_active_cells;
	return active_cells;
}


int move_nodes(CELLS * grid, int * active_cells, int num_active_cells,
 int * nc, double * l, MAT * nodes)
{

	CELL ** grid_cells = grid->cells;

	// Only loop over the active_cells
	int i = 0;
	int j = 0;

	int i_new = 0;
	int j_new = 0;


	int k = 0;
	int k_new = 0;
	int cell_number = 0;
	int dim = nodes->n;

	NODELIST * current_p = NULL;
	NODELIST * previous_p = NULL;

	double x,y,z;

	if ( dim == 2)
	{

		int nx = nc[0];
		int ny = nc[1];

		for ( int i = 0 ; i < num_active_cells ; i++)
		{

			cell_number = active_cells[i];

			// translate into each cell number into its i,j grid reference

			i = (int) (floor(cell_number/nx));
			j = (int) (cell_number - i*(nx));


			current_p = grid_cells[i][j].nodes;

			if ( current_p == NULL)
			{
				// remove cell from active list 
			}
			previous_p = current_p;

			while ( current_p != NULL)
			{

				// get coordinates 
				x = nodes->me[current_p->node_number][0];
				y = nodes->me[current_p->node_number][1];

				// Find grid reference
				i_new = (int)floor(x*nc[0]/l[0]);
				j_new = (int)floor(y*nc[1]/l[1]);

				// check if correct grid reference
				if (( i_new != i ) ||  (j_new != j))
				{
					// remove point from current grid nodelist
					deleteList(previous_p);
					// add point to current node list 
					insertList(&grid_cells[i_new][j_new],current_p);
					// set current point to the previous point
					current_p = previous_p;


				}else{
					// update previous point to current point
					previous_p = current_p;
				}

				// move on to next point in list
				current_p = current_p->next;


			}






		}






	}



	return num_active_cells;
}



int neighbour_RangeSearch(IVEC * neighbours, int cell,
int dim,  double * x, double * range, RANGE_TYPE range_type)
{
	// Find all nodes that are within this range 
	int num_neighbours = 0;


	// Find the possible cells


	switch ( range_type)
	{
		case(RADIAL_SEARCH):
		{
			// iterate over each cell


			// Find each point within the range of point x
			// as measured by the euclidean norm




			break;
		}
		case(RECTANGULAR_SEARCH):
		{

			// iterate over each cell

			break;
		}

	}



	// check which possible cells it could be 

	return num_neighbours;
}