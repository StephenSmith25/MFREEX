#include "cellSearch.h"


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
