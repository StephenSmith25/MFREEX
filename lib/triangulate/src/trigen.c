/*!
 This function takes in a planar-straight length graph defined in the .nodes and .segs files
 and returns the triangulion, subject to the options given in char[] options
*/

#include <trigen.h>
#include <connect_segments.h>

void report(io, markers, reporttriangles, reportneighbors, reportsegments,
		reportedges, reportnorms)
	struct triangulateio *io;
	int markers;
	int reporttriangles;
	int reportneighbors;
	int reportsegments;
	int reportedges;
	int reportnorms;
{
	int i, j;

	for (i = 0; i < io->numberofpoints; i++) {
		printf("Point %4d:", i);
		for (j = 0; j < 2; j++) {
			printf("  %.6g", io->pointlist[i * 2 + j]);
		}
		if (io->numberofpointattributes > 0) {
			printf("   attributes");
		}
		for (j = 0; j < io->numberofpointattributes; j++) {
			printf("  %.6g",
					io->pointattributelist[i * io->numberofpointattributes + j]);
		}
		printf("\n");
	}


	if (reporttriangles || reportneighbors) {
		for (i = 0; i < io->numberoftriangles; i++) {
			if (reporttriangles) {
				printf("Triangle %4d points:", i);
				for (j = 0; j < io->numberofcorners; j++) {
					printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
				}
				if (io->numberoftriangleattributes > 0) {
					printf("   attributes");
				}
				for (j = 0; j < io->numberoftriangleattributes; j++) {
					printf("  %.6g", io->triangleattributelist[i *
							io->numberoftriangleattributes + j]);
				}
				printf("\n");
			}
			if (reportneighbors) {
				printf("Triangle %4d neighbors:", i);
				for (j = 0; j < 3; j++) {
					printf("  %4d", io->neighborlist[i * 3 + j]);
				}
				printf("\n");
			}
		}
		printf("\n");
	}

	if (reportsegments) {
		for (i = 0; i < io->numberofsegments; i++) {
			printf("Segment %4d points:", i);
			for (j = 0; j < 2; j++) {
				printf("  %4d", io->segmentlist[i * 2 + j]);
			}
			if (markers) {
				printf("   marker %d\n", io->segmentmarkerlist[i]);
			} else {
				printf("\n");
			}
		}
		printf("\n");
	}

	if (reportedges) {
		for (i = 0; i < io->numberofedges; i++) {
			printf("Edge %4d points:", i);
			for (j = 0; j < 2; j++) {
				printf("  %4d", io->edgelist[i * 2 + j]);
			}
			if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
				for (j = 0; j < 2; j++) {
					printf("  %.6g", io->normlist[i * 2 + j]);
				}
			}
			if (markers) {
				printf("   marker %d\n", io->edgemarkerlist[i]);
			} else {
				printf("\n");
			}
		}
		printf("\n");
	}
}

TRIANGLE  *  trigen(char * options, char * fileName)
{

	char nodesFile[20] ;
	char segsFile[20] ;


	strcpy(nodesFile,fileName);
	strcpy(segsFile,fileName);
	strcat(nodesFile,".nodes");
	strcat(segsFile,".segs");

	printf("%s\n",nodesFile );
	printf("%s\n",segsFile );
	int i;


	struct triangulateio in, out;


	

	/*  Input structure */
	char buf[64];
	char s[2] = " ";

	char * token;

	FILE * fp = fopen(nodesFile,"r");
	if (fgets(buf,sizeof(buf),fp) != NULL)
		token = strtok(buf,s);
	in.numberofpoints = atoi(token);
	token = strtok(NULL,s);
	int isAttribute = atoi(token);


	printf("is attruibute = %d \n", isAttribute);
	if ( isAttribute == 1){
		in.numberofpointattributes = 1;
		in.pointattributelist = (double* )  malloc(in.numberofpoints*sizeof(double));
	}else{
		in.numberofpointattributes = 0;
		in.pointattributelist = NULL;
	}
	in.pointlist = malloc(in.numberofpoints*2*sizeof(double));
	for ( i = 0 ; i < in.numberofpoints ;i++){
		if (fgets(buf,sizeof(buf),fp) != NULL)
			token = strtok(buf,s);
		in.pointlist[2*i] = atof(token);
		token = strtok(NULL,s);
		in.pointlist[2*i+1] = atof(token);

		if ( isAttribute == 1){
			token = strtok(NULL,s);
			in.pointattributelist[i] = atof(token);
		}

	}
	fclose(fp);

	fp = fopen(segsFile,"r");
	if (fgets(buf,sizeof(buf),fp) != NULL)
		in.numberofsegments = atoi(buf);

	printf("number of segments = %d \n",in.numberofsegments);
	in.segmentlist = malloc(2*in.numberofsegments*sizeof(int));
	in.segmentmarkerlist = malloc(in.numberofsegments*sizeof(int));
	in.pointmarkerlist = (int * ) NULL;


	printf("got here \n");
	for ( i = 0 ; i < in.numberofsegments ;i++){
		if (fgets(buf,sizeof(buf),fp) != NULL)
			token = strtok(buf,s);
		in.segmentlist[2*i] = atoi(token);
		token = strtok(NULL,s);
		in.segmentlist[2*i+1] = atoi(token);
		token = strtok(NULL,s);
		if ( token != NULL)
			in.segmentmarkerlist[i] = atoi(token);

	}


	in.numberofregions = 0;
	in.regionlist = (REAL *) NULL;
	//holes
	in.numberofholes = 0;
	in.holelist = (double *) NULL;

	//report(&in, 1, 0, 0, 1, 0, 0);

	/* Make necessary initializations so that Triangle can return a */
	/*   triangulation in `out' and a voronoi diagram in `vorout'.  */

	out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
	/* Not needed if -N switch used or number of point attributes is zero: */
	out.pointattributelist = (double *) NULL;
	out.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
	out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
	/* Not needed if -E switch used or number of triangle attributes is zero: */
	out.triangleattributelist = (double *) NULL;
	out.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
	/* Needed only if segments are output (-p or -c) and -P not used: */
	out.segmentlist = (int *) NULL;
	/* Needed only if segments are output (-p or -c) and -P and -B not used: */
	out.segmentmarkerlist = (int *) NULL;


	
	/* Triangulate the points.  Switches are chosen to read and write a  */
	/*   PSLG (p), preserve the convex hull (c), number everything from  */
	/*   zero (z), assign a regional attribute to each element (A), and  */
	/*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
	/*   neighbor list (n).                                              */


	/*  Call triangulate */
	triangulate(options,&in,&out, (struct trianglulateio * ) NULL);

	report(&out, 1, 0, 0, 0, 0, 0);



	// find boundary
	int * connected_segments = NULL;
	connect_segments(out.segmentlist, &connected_segments, out.numberofsegments);


	// prepare output points
	int * boundaryPoints = malloc((out.numberofsegments)*sizeof(int));

	for (i = 0; i < out.numberofsegments; i++) {
		boundaryPoints[i] = connected_segments[2*i] -1 ;
		/* code */
	}

	TRIANGLE * triangular_mesh = malloc(1*sizeof(TRIANGLE));

	if ( isAttribute == 1)
	{
		triangular_mesh->temperatures = out.pointattributelist;

	}


	triangular_mesh->points = out.pointlist;
	triangular_mesh->boundary = boundaryPoints;
	triangular_mesh->num_points = out.numberofpoints;
	triangular_mesh->num_boundary_points = out.numberofsegments;
	triangular_mesh->pointmarkers = out.pointmarkerlist;
	triangular_mesh->triangles = out.trianglelist;
	triangular_mesh->num_triangles = out.numberoftriangles;




	// end of program, free memory !

	// free(connected_segments);
	// free(in.pointlist);
 // 	free(in.pointattributelist);
 //  	free(in.segmentmarkerlist);
 //  	free(in.segmentlist);
 //  	free(in.regionlist);


  	return triangular_mesh;
	



}


	//in.pointlist = (REAL *) NULL;
	//in.pointattributelist = (REAL * ) NULL;
	//in.pointmarkerlist = (int * ) NULL;
	// in.trianglelist = (int * ) NULL;
	// in.triangleattributelist = (REAL*) NULL;
	// in.trianglearealist = (REAL *) NULL;
	// in.neighborlist = (int *) NULL;
	// in.segmentlist = (int *) NULL;
	// in.segmentmarkerlist = (int * ) NULL;
	// in.holelist = (REAL*) NULL;
	// in.regionlist = (REAL*) NULL;
	// in.edgelist = (int *) NULL;
	// in.edgemarkerlist = (int *) NULL;
	// in.normlist = (REAL * ) NULL;