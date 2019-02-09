
#include "input/read_input_points.h"

const int DIM = 2;

int read_input_points(char * fileName, double ** points_out, int ** boundary_out, 
	int ** boundary_markers_out, double ** attributes_out, int* num_nodes_out, int *num_boundary_out)
{

	FILE * fp;

	// deliminator 
	char delim[2] = " ";
	char * token;
	// boundary 
	char boundaryFile[50];
	char buffer[50];
	strcpy(boundaryFile,fileName);
	strcat(boundaryFile,".boundary");

	fp = fopen(boundaryFile,"r");
	fgets(buffer, sizeof(buffer),fp);
	token = strtok(buffer,delim );
	int num_boundary = atoi(token);



	
	int * boundary = malloc(num_boundary * sizeof(int));
	int * boundary_markers = malloc(num_boundary * sizeof(int));


	for ( int i = 0 ; i < num_boundary ; i++)
	{
		fgets(buffer, sizeof(buffer),fp);
		token = strtok(buffer,delim );
		boundary[i] = atoi(token)-1;
		token = strtok(NULL, delim);
		boundary_markers[i] = atoi(token);

	}

	*boundary_out = boundary;
	*boundary_markers_out = boundary_markers;
	*num_boundary_out = num_boundary;
	// nodes 

	fclose(fp);

	char nodesFile[50];
	strcpy(nodesFile,fileName);
	strcat(nodesFile,".nodes");

	fp = fopen(nodesFile,"r");
	fgets(buffer, sizeof(buffer),fp);
	token = strtok(buffer,delim );
	int num_nodes = atoi(token);
	token = strtok(NULL, delim);

	int is_attribute = atoi(token);	

	double * attributes = NULL;
	if ( is_attribute == 1)
	{
		attributes = (double*) malloc(num_nodes * sizeof(double));
	}



	double * points = malloc(DIM * num_nodes *  sizeof(double) );

	for ( int i = 0 ; i < num_nodes ; i++)
	{
		fgets(buffer, sizeof(buffer),fp);
		token = strtok(buffer,delim);

		points[2*i] = atof(token);
		token = strtok(NULL,delim);
		points[2*i+1] = atof(token);

		if ( is_attribute)
		{
			token = strtok(NULL,delim);

			attributes[i] = atof(token);
		}

	}

	fclose(fp);

	*num_nodes_out = num_nodes;
	*attributes_out = attributes;
	*points_out = points;





	return 0;
}