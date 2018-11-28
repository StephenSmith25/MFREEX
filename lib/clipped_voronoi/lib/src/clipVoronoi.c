#include "clipVoronoi.h"



int clipVoronoi(gpc_polygon*** voronoi, double points[], int boundary[], int numPoints, int numBoundary)
{

  FILE * fp;

  // create the clipping polygon
  gpc_vertex clipVertex[numBoundary];


  for (int i = 0; i < numBoundary; i++) {
    clipVertex[i].x = points[2*boundary[i]];
    clipVertex[i].y = points[2*boundary[i]+1];
  }

  gpc_vertex_list clip_list = {.num_vertices = numBoundary, .vertex = clipVertex};
  gpc_polygon clip_polygon = {.num_contours=1, .hole = NULL, .contour = &clip_list};

  fp = fopen("cells.txt","w");
  gpc_write_polygon(fp, 0, &clip_polygon);
  fclose(fp);
  printf("CLIPPING VORONOI\n");

  // clip boundary cells
  for (size_t i = 0; i < numPoints; i++) {

    int indx = i;
    gpc_polygon * clipped_polygon = malloc(1*sizeof(gpc_polygon));
    gpc_polygon_clip(GPC_INT,&clip_polygon,(*voronoi)[indx],clipped_polygon);
    gpc_free_polygon((*voronoi)[indx]);
    free((*voronoi)[indx]);
    (*voronoi)[indx] = clipped_polygon;
    /* code */
  }

  printf("FINISHED CLIPPING VORONOI\n");

 //write cells to file

  // char fileName[256];
  // for (int i = 0; i < numPoints; i++) {
  //   // write cell to file
  //   	snprintf(fileName,sizeof(fileName),"Cells/%d.txt",i);
  //   	// write clipped cell to files
  //   	fp = fopen(fileName,"w");
  //   	if ( fp != NULL)
  //   		gpc_write_polygon(fp, 0, (*voronoi)[i]);
  //   	fclose(fp);
  
  // }
  
  


  // put voronoi into a vertex point list and
  //create an array to store how many verticies each polygon has

  
  return 0;
}
