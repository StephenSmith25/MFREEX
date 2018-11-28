
#include "callVoronoi.h"



#define JCV_ATAN2 atan2
#define JCV_FLT_MAX 1.7976931348623157E+308

gpc_polygon ** callVoronoi(double points[], int num_points){

  // put points into jcv points structure and then call voronoi


  jcv_point voronoi_points[num_points] ;

  for (int i = 0; i < num_points; i++) {
    voronoi_points[i].x = points[2*i];
    voronoi_points[i].y = points[2*i+1];
  }


  // allocate memory for voronoi diagram
  jcv_diagram diagram;
  const jcv_site * sites = NULL;
  jcv_graphedge * graph_edge;
  memset(&diagram,0,sizeof(jcv_diagram));

  printf("GENERAING VORONOI DIAGRAM \n");
  jcv_diagram_generate(num_points,(const jcv_point *)voronoi_points,NULL,&diagram);
  printf("FINISHED GENERATING VORONOI DIAGRAM \n");

  // get voronoi sites
  sites = jcv_diagram_get_sites(&diagram);

  // move voronoi into gpc_polygon structure, where each polygon is stored
  // as a group of verticies
  gpc_polygon ** voronoi = malloc(diagram.numsites*sizeof(gpc_polygon*));


  for (size_t i = 0; i < diagram.numsites; i++) {
    voronoi[i] = malloc(1*sizeof(gpc_polygon));
    voronoi[i]->hole = NULL;
    voronoi[i]->num_contours = 1;
    voronoi[i]->contour = malloc(1*sizeof(gpc_vertex_list));
    voronoi[i]->contour->num_vertices = 0;
    voronoi[i]->contour->vertex = NULL;
  }

  // loop over sites and set up output polygons
  for (size_t i = 0; i < diagram.numsites; i++) {

    // number of verticies in each cell
    int numVert = 0;

    graph_edge = sites[i].edges;
    int indx = sites[i].index;


    // loop over cell to get number of verticies
    while(graph_edge){
      ++numVert;
      graph_edge = graph_edge->next;
    }

    // assign memory for the polygons
    voronoi[indx]->contour->num_vertices = numVert;
    voronoi[indx]->contour->vertex = malloc(numVert*sizeof(gpc_vertex));

    // loop over edges
    graph_edge = sites[i].edges;

    for ( int k = 0 ; k < numVert ; k++){
      voronoi[indx]->contour->vertex[k].x = graph_edge->pos[1].x;
      voronoi[indx]->contour->vertex[k].y = graph_edge->pos[1].y;
      graph_edge = graph_edge->next;
    }



  }


  jcv_diagram_free(&diagram);

  return voronoi;



  return 0;


}
