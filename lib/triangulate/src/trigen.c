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

struct triangulateio *  trigen(double ** output_points, int ** boundary, char * options, char * fileName, int * numPoints, int * numBoundary, int ** pointmarkers, double ** temperatures ){

	char nodesFile[20] ;
	char segsFile[20] ;


	strcpy(nodesFile,fileName);
	strcpy(segsFile,fileName);
	strcat(nodesFile,".nodes");
	strcat(segsFile,".segs");

	printf("%s\n",nodesFile );
	printf("%s\n",segsFile );
	int i;


	tri in = malloc(sizeof(struct triangulateio ));
	tri out = malloc(sizeof(struct triangulateio));
		/*  Input structure */
	char buf[64];
	char s[2] = " ";

	char * token;

	FILE * fp = fopen(nodesFile,"r");
	if (fgets(buf,sizeof(buf),fp) != NULL)
		token = strtok(buf,s);
	in->numberofpoints = atoi(token);
	token = strtok(NULL,s);
	int isAttribute = atoi(token);


	printf("is attruibute = %d \n", isAttribute);
	if ( isAttribute == 1){
		in->numberofpointattributes = 1;
		in->pointattributelist = (double* )  malloc(in->numberofpoints*sizeof(double));
	}else{
		in->numberofpointattributes = 0;
		in->pointattributelist = NULL;
	}
	in->pointlist = malloc(in->numberofpoints*2*sizeof(double));
	for ( i = 0 ; i < in->numberofpoints ;i++){
		if (fgets(buf,sizeof(buf),fp) != NULL)
			token = strtok(buf,s);
		in->pointlist[2*i] = atof(token);
		token = strtok(NULL,s);
		in->pointlist[2*i+1] = atof(token);

		if ( isAttribute == 1){
			token = strtok(NULL,s);
			in->pointattributelist[i] = atof(token);
		}

	}
	fclose(fp);
	fp = fopen(segsFile,"r");
	if (fgets(buf,sizeof(buf),fp) != NULL)
		in->numberofsegments = atoi(buf);
	in->segmentlist = malloc(2*in->numberofsegments*sizeof(int));
	in->segmentmarkerlist = malloc(in->numberofsegments*sizeof(int));
	for ( i = 0 ; i < in->numberofsegments ;i++){
		if (fgets(buf,sizeof(buf),fp) != NULL)
			token = strtok(buf,s);
		in->segmentlist[2*i] = atoi(token);
		token = strtok(NULL,s);
		in->segmentlist[2*i+1] = atoi(token);
		token = strtok(NULL,s);
		if ( token != NULL)
			in->segmentmarkerlist[i] = atoi(token);

	}
	in->numberofregions = 0;
	in->regionlist = (double *) NULL;
		//holes
	in->numberofholes = 0;
	in->holelist = (double *) NULL;

	// set of options inserted
	//printf("Input point set:\n");
	//report(in, 1, 0, 0, 1, 0, 0);

	/* Make necessary initializations so that Triangle can return a */
	/*   triangulation in `out' and a voronoi diagram in `vorout'.  */

	out->pointlist = (double *) NULL;            /* Not needed if -N switch used. */
	/* Not needed if -N switch used or number of point attributes is zero: */
	out->pointattributelist = (double *) NULL;
	out->pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
	out->trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
	/* Not needed if -E switch used or number of triangle attributes is zero: */
	out->triangleattributelist = (double *) NULL;
	out->neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
	/* Needed only if segments are output (-p or -c) and -P not used: */
	out->segmentlist = (int *) NULL;
	/* Needed only if segments are output (-p or -c) and -P and -B not used: */
	out->segmentmarkerlist = (int *) NULL;
	out->edgelist = (int *) NULL;             /* Needed only if -e switch used. */
	out->edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */
	out->holelist = (double *) NULL;
	out->regionlist = (double * ) NULL;
	/* Triangulate the points.  Switches are chosen to read and write a  */
	/*   PSLG (p), preserve the convex hull (c), number everything from  */
	/*   zero (z), assign a regional attribute to each element (A), and  */
	/*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
	/*   neighbor list (n).                                              */

	/*  Failed at triangulate */
	triangulate(options,in,out, (struct trianglulateio *) NULL);

	//printf("passed triangulate \n");
	//report(out, 1, 0, 0, 0, 0, 0);



	// find boundary
	int * connected_segments = NULL;
	connect_segments(out->segmentlist, &connected_segments, out->numberofsegments);


	// prepare output points
	int * boundaryPoints = malloc((out->numberofsegments)*sizeof(int));

	for (i = 0; i < out->numberofsegments; i++) {
		boundaryPoints[i] = connected_segments[2*i] -1 ;
		/* code */
	}


	*output_points = out->pointlist;
	*boundary = boundaryPoints;
	*numPoints = out->numberofpoints;
	*numBoundary = out->numberofsegments ;
	*pointmarkers = out->pointmarkerlist;

	if ( isAttribute == 1)
	{
		*temperatures = out->pointattributelist;

	}

	(*pointmarkers) = malloc(out->numberofpoints*sizeof(int));

	for ( i = 0 ; i < out->numberofpoints ; i++)
	{
		(*pointmarkers)[i] = out->pointmarkerlist[i];

	}

	// end of program, free memory !

	free(connected_segments);
	free(in->pointlist);
 	free(in->pointattributelist);
  	free(in->segmentmarkerlist);
  	free(in->segmentlist);
  	free(in->regionlist);
  	free(in);


  	return out;
	



}
