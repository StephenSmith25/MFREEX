#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "matrix2.h"
#include "meshgrid.h"
#include "generate_voronoi.h"
#include "meshgrid.h"
#include <sys/time.h>
#include "Material/material.h"
#include "Integration/SCNI/generate_scni.h"
#include "Integration/SCNI/generate_mscni.h"
#include "mls_shapefunction.h"
#include "setDomain.h"
#include "smoothstep.h"
#include "Force/Internal/internalForce_hyperelastic.h"
#include "Boundary/Displacement/essential_boundary.h"
#include "mat2csv.h"
#include "trigen.h"
#include "Boundary/getBoundary.h"
#include "Boundary/iv_addNode.h"
#include "Boundary/Traction/Cavity/cavityVolume.h"
#include "Boundary/Traction/Cavity/flowRate.h"
#include "Boundary/Traction/Pressure/pressure_load.h"
#include <math.h>
#include "Deformation/poldec.h"
#include "Boundary/Contact/contactDetection.h"
#include "Boundary/Displacement/setUpBC.h"
#include "Boundary/Displacement/enforceBC.h"
#include "Integration/SCNI/scni_update_B.h"
#include "Integration/material_point.h"
#include "Integration/mass_matrix.h"
#include "Integration/DomainMaterialPoint.h"
#include "cellSearch.h"
#include <stdbool.h>
#include <pthread.h>
#include "internal_force_mooney.h"
#include "input/read_input_mesh.h"
#include "input/read_config.h"
#include "Boundary/Displacement/DisplacementBC.h"



#define NUM_SEARCHES 1e6


int nnx = 100;
int nny = 100;




int main(int argc, char** argv) {



	
	/*////////////////////////////////////////////////////////// */
	/*                                                           */
	/*			           Create Geometry                       */
	/*                                                           */
	/*////////////////////////////////////////////////////////// */

  

	double X_MAX = 10;
	double X_MIN = 0;
	double Y_MAX = 10;
	double Y_MIN = 0;


	double testpoint[2] = {5,5};

	int dim = 2;
	int * boundaryNodes ;
	int num_boundary;


	// CREATE NODES
	MAT * xI = meshgrid(X_MIN,X_MAX,nnx,
		Y_MIN,Y_MAX,nny, &boundaryNodes, &num_boundary);

	FILE * fp;

	int numnodes = xI->m;
	fp = fopen("nodes.csv","w");
	for ( int i =0 ; i < numnodes ; i++)
	{
		for ( int k = 0 ; k < dim ; k++)
		{
			fprintf(fp, "%lf",xI->me[i][k]);

			if ( k < dim - 1)
			{
				fprintf(fp,",");
			}

		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	

	


	// CREATE NODES
	BOUNDING_BOX * bounding_box = create_bounding_box(X_MIN, X_MAX+1e-9,
	Y_MIN,Y_MAX+1e-9,0, 0);

	double cell_size[2] = {0.4,0.4};

	MAT * xI_copy = m_copy(xI,MNULL);
	CELLS * cells = create_cells(bounding_box, cell_size, dim, xI_copy);

	fp = fopen("search_cells.csv","w");
	for ( int i = 0 ; i < cells->ny ; i++)
	{

		for ( int j = 0 ; j < cells->nx ; j++)
		{

		fprintf(fp,"%lf,%lf,%lf,%lf\n",cells->cells[i][j].x[0],cells->cells[i][j].x[1],
			cells->cells[i][j].y[0],cells->cells[i][j].y[1]);
		
		}

	}




	fclose(fp);

	int num_active_cells = 0;
	active_cell * active_cells = get_active_cells(cells, &num_active_cells);
	write_active_cells("active_cells.csv",active_cells);
	
	struct timeval start, end;


	// NAIVE SEARCH
	double r = 0.5;



	int num_neighbours = 0;
	int * neighbours = malloc(20*sizeof(int));
	gettimeofday(&start, NULL);

	for (int k = 0 ; k < NUM_SEARCHES ; k++)
	{
		num_neighbours = 0;


		for ( int i = 0 ; i < numnodes ; i++)
		{

			double x = xI->me[i][0];
			double y = xI->me[i][1];

			double r_temp = sqrt(pow(x -testpoint[0],2) + pow(y-testpoint[1],2));

			if ( r_temp < r)
			{
				neighbours[num_neighbours] = i;
				++num_neighbours;
			}

		}
	}
	gettimeofday(&end, NULL);
	double time_taken =((end.tv_sec  - start.tv_sec) * 1000000u + 
         		end.tv_usec - start.tv_usec) / 1.e6;

	printf("time taken with naive= %lf \n", time_taken);



	// create the material point
	MATERIAL_POINT * MP= malloc(1*sizeof(MATERIAL_POINT));
	//MP->MI = m_get(2,2);
	MP->invMI = m_get(2,2);
	MP->invMI->me[0][0] = 1.00/(r*r);
	MP->invMI->me[1][1] = 1.00/(r*r);
	MP->neighbours = iv_get(80);
	MP->coords_n_1 = malloc(2*sizeof(double));
	MP->coords_n_1[0] = testpoint[0];
	MP->coords_n_1[1] = testpoint[1];


	gettimeofday(&start, NULL);


	for (int k = 0 ; k < NUM_SEARCHES ; k++)
	{

		RangeSearchMaterialPoint(MP, xI,cells);

	}

	gettimeofday(&end, NULL);
	double time_taken_1 =((end.tv_sec  - start.tv_sec) * 1000000u + 
         		end.tv_usec - start.tv_usec) / 1.e6;

	// printf("APPROACH 1\n ---------\n");

	// for ( int k = 0 ; k < num_neighbours ; k++)
	// {
	// 	printf("%d,",neighbours[k]);

	// }

	// printf("\n");


	// printf("APPROACH 2\n ---------\n");

	// for ( int k = 0 ; k < num_neighbours ; k++)
	// {
	// 	printf("%d,",MP->neighbours->ive[k]);

	// }
	// printf("\n");

	printf("time taken with cells = %lf \n", time_taken_1);



	printf("number of neibhours = %d \n",num_neighbours);
	printf("number_of_neighbours = %d \n",MP->num_neighbours);


	exit(0);
}
