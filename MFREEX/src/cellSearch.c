#include "cellSearch.h"

static inline double sq_distance(double *point_1, double * point_2, int  dim)
{
	double distance = 0;
	for (int i = 0; i < dim; ++i)
	{
		distance += pow(point_2[i] - point_1[i], 2);
		/* code */
	}

	return sqrt(distance);

}



static inline void deleteList(NODELIST ** nodelist, NODELIST * delete, NODELIST * prev)
{

	if (*nodelist == delete){
		*nodelist = delete->next;

	}else if ( delete->next != NULL){
		prev->next = delete->next;
	}
	else{
		prev->next = NULL;
	}

}
static inline void insertList(CELL * cell, NODELIST * newpoint)
{
	// link the new item to point ot the head of the list
	newpoint->next = (cell)->nodes;

	// set of head of the list be our new item
	(cell)->nodes = newpoint;



}


CELLS * create_cells(BOUNDING_BOX * bounding_box, double * CELL_SIZE, int dim, MAT * NODES)
{



	CELLS * cells = malloc(1*sizeof(CELLS));


	if ( dim == 2)
	{


		// GET CELL DIMENSIONS
		double CX = CELL_SIZE[0];
		double CY = CELL_SIZE[1];
		cells->CX = CX;
		cells->CY = CY;



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
		cells->cells = malloc(ny*sizeof(CELL*));
		for ( int i =0 ; i < ny ; i++)
		{
			cells->cells[i] = malloc(nx*sizeof(CELL));
		}

		for ( int j = 0 ; j < ny ; j++ )
		{
			for (int i = 0 ; i < nx ; i++)
			{
				cells->cells[j][i].x = malloc(2*sizeof(double));
				cells->cells[j][i].y = malloc(2*sizeof(double));


				/*Create cell limits*/
				cells->cells[j][i].x[0] = i*CX;
				cells->cells[j][i].x[1] = (i+1)*CX;

				cells->cells[j][i].y[0] = j*CY;
				cells->cells[j][i].y[1] = (j+1)*CY;

				// CREATE CELL NODE LIST
				cells->cells[j][i].nodes = NULL;

				cells->cells[j][i].num_nodes = 0;



				}


			}

			for ( int k = 0 ; k < NODES->m ; k++)
			{
				
					double x = NODES->me[k][0];
					double y = NODES->me[k][1];


					int j = (int)floor(x/CX);
					int i = (int)floor(y/CY);


					/*Create nodelist object */
					NODELIST * node_i = (NODELIST * )malloc(1*sizeof(NODELIST));
					node_i->node_number = k;

					insertList(&cells->cells[i][j], node_i);	
					++cells->cells[i][j].num_nodes ;



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

static inline void deleteCell(active_cell ** active_cells, int cell_to_delete)
{
	
	// check if active_cells is already in this list
	active_cell * previous = *active_cells;
	active_cell * current = *active_cells;
	active_cell * TEMP = *active_cells;

	int cell_number = 0;



	while ( current != NULL) 
	{

		cell_number = current->cell_number;
		if ( cell_number == cell_to_delete)
		{
			// check if at head
			if ( current == previous)
			{
				TEMP = current->next;
				free(current);
				*active_cells = TEMP;

			}

			// check if at tail
			else if( current->next == NULL)
			{
				free(current);
				previous->next = NULL;
			}


			// else normal
			else
			{
				TEMP = current->next;
				free(current);
				previous->next = TEMP;

			}

			return;



		}else{
			previous = current;
			current = current->next;
		}

	}
	return;

}

static inline void addCell(active_cell ** active_cells , int cell_to_add )
{

	// check if active_cells is already in this list
	active_cell * previous = *active_cells;
	active_cell * current = *active_cells;

	int cell_number = 0;
	int stop = 0;
	while (( current != NULL) && (stop != 1))
	{

		cell_number = current->cell_number;


		if ( cell_to_add == cell_number) {
			stop = 1;
			return;
		}else if ( cell_number > cell_to_add ) {

			/* Create the new cell */
			active_cell * active_cell_i = malloc(1*sizeof(active_cell));
			active_cell_i->cell_number = cell_to_add;

			// check if head
			if ( current == previous)
			{

				//active_cell_i->next = current;
				*active_cells = active_cell_i;
				active_cell_i->next = current;

			// else add somewhere inbetween the head and tail 
			}else{
				previous->next = active_cell_i;
				active_cell_i->next = current;
			}
		
			stop = 1;


		}else if ( current->next == NULL) {
			/* Create the new cell */
			active_cell * active_cell_i = malloc(1*sizeof(active_cell));
			active_cell_i->cell_number = cell_to_add;

			// add to tail
			active_cell_i->next = NULL;
			current->next = active_cell_i;

			stop = 1;

		}else{
			previous = current;
			current = current->next;

		}

	
	}
	return;

}

int write_active_cells(char * filename, active_cell * active_cells)
{
	FILE * fp = fopen(filename,"w");
	active_cell * activeCell = active_cells;

	while ( activeCell != NULL)
	{	
		fprintf(fp,"%d \n",activeCell->cell_number);
		activeCell = activeCell->next;
	}
	fclose(fp);

	return 0;
}
active_cell * get_active_cells(CELLS * cells, int * num_active_cells)
{

	int i,j,k; 

	int count_num_active_cells = 0;
	int nx = cells->nx;
	int ny = cells->ny;

	active_cell * active_cells = NULL;
	active_cell * previous = NULL;

	for ( int i = 0 ; i < cells->ny ; i++)
	{

		for ( int j = 0 ; j < cells->nx ; j++)
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



int move_nodes(CELLS * grid, active_cell ** active_cells, double * l, MAT * nodes)
{

	CELL ** grid_cells = grid->cells;

	// Only loop over the active_cells
	int i = 0;
	int j = 0;

	int i_new = 0;
	int j_new = 0;

	int nx = grid->nx;
	int ny = grid->ny;

	int k = 0;
	int k_new = 0;
	int cell_number = 0;
	int cell_number_new = 0;
	int dim = nodes->n;

	NODELIST * current_p = NULL;
	NODELIST * previous_p = NULL;

	double x,y,z;
	int count = 0;

	if ( dim == 2)
	{

		active_cell * active_cell_i = *active_cells;
		while ((active_cell_i != NULL ) && ( count < 5000))
		{
			++count;

			cell_number = active_cell_i->cell_number;
			// translate into each cell number into its i,j grid reference

			i = (int) (floor(cell_number/nx));
			j = (int) (cell_number - i*(nx));



			current_p = grid_cells[i][j].nodes;

			/* If node list is empty delete cell from active cells*/
			if ( current_p == NULL)
			{

				deleteCell(active_cells, cell_number);


			}
			previous_p = current_p;


			while ( current_p != NULL)
			{

				// get coordinates 
				x = nodes->me[current_p->node_number][0];
				y = nodes->me[current_p->node_number][1];


				// Find grid reference
				j_new = (int)floor(x/l[0]);
				i_new = (int)floor(y/l[1]);


				// check if correct grid reference
				if (( i_new != i ) ||  (j_new != j))
				{

					// Add new cell to active list
					// if already a member, cell wont be added 
					cell_number_new = (i_new)*nx + j_new;
					addCell(active_cells, cell_number_new);
					// remove point from current grid nodelist
					deleteList(&grid_cells[i][j].nodes,current_p, previous_p);
					// add point to new node list 
					insertList(&grid_cells[i_new][j_new],current_p);


					// set current point to the previous point
				//	current_p = previous_p->next;	

					if ( current_p == previous_p)
					{
						//grid_cells[i][j].nodes = previous_p->next;
						current_p = grid_cells[i][j].nodes;
						previous_p = current_p;
					}else{
						current_p = previous_p->next;	
	
					}

				}else{

					// update previous point to current point
					previous_p = current_p;
					// move on to next point in list
					current_p = current_p->next;


				}



			}




		active_cell_i = active_cell_i->next;

		}






	}



	return 0;
}



int neighbour_RangeSearch(IVEC * neighbours, CELLS * cells,
double * x, double range,  MAT * nodes)
{
	// NOTE - IMPLEMENTED JUST FOR 2D AT THE MOMENT - END NOTE 

	// Find all nodes that are within this range 
	int num_neighbours = 0;
	int dim = nodes->n;
	
	// convert x into its i,j components
	int j = (int)floor(x[0]/cells->CX);
	int i = (int)floor(x[1]/cells->CY);

	// number of clels in each direction
	int nx = cells->nx;
	int ny = cells->ny;



	NODELIST * current_p = NULL;
	double * node_check = NULL;
	double distance = 0;
	for ( int indx_j = j-1 ; indx_j <= j+1 ; indx_j++)
		for ( int indx_i = i-1 ; indx_i <= i+1 ; indx_i++)
	{


		if (( indx_j >= 0) && (indx_j < nx) && ( indx_i >= 0) && (indx_i < ny))
		{
			current_p = cells->cells[indx_i][indx_j].nodes;
			while ( current_p != NULL)
			{

				node_check = nodes->me[current_p->node_number];

				distance = sq_distance(x,node_check, dim);

				if ( distance <= range)
				{
					neighbours->ive[num_neighbours] = current_p->node_number;
					++num_neighbours;
				}


				current_p = current_p->next;

			}


		}else{
			// skip
		}


	}

	return num_neighbours;
}