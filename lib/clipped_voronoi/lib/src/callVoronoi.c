
#include "callVoronoi.h"
#define JC_VORONOI_IMPLEMENTATION
#define JCV_REAL_TYPE double
#define JCV_ATAN2 atan2
#define JCV_FLT_MAX 1.7976931348623157E+308
#define JCV_CEIL ceil
#define JCV_FLOOR floor
#define JCV_FABS fabs
#include "jc_voronoi.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

static inline int is_equal_double(double * a, double * b, int dim){
  double tol = 1e-14;


  if (     (fabs(a[0] - b[0]) < tol) && ( fabs(a[1] - b[1]) < tol ) )
  {
    return 1;
  }

  return 0;

}

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

  jcv_diagram_generate(num_points,(const jcv_point *)voronoi_points,NULL,&diagram);



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

  //diagram.numsites
  // loop over sites and set up output polygons
  for (int i = 0; i < diagram.numsites ; i++) {


    // number of verticies in each cell
    int numVert = 0;

    graph_edge = sites[i].edges;
    int indx = sites[i].index;


    // loop over cell to get number of verticies

    double v_1[2] = {-1000,1000};
    double v_2[2] = {-1000,-1000};


    while(graph_edge){
      v_2[0] = graph_edge->pos[1].x;
      v_2[1] = graph_edge->pos[1].y;

      if ( is_equal_double(v_1, v_2, 2) != 1){
        ++numVert;
        v_1[0] = v_2[0];
        v_1[1] = v_2[1];
      }
      graph_edge = graph_edge->next;

    }


    // assign memory for the polygons
    voronoi[indx]->contour->num_vertices = numVert;
    voronoi[indx]->contour->vertex = malloc(numVert*sizeof(gpc_vertex));

    // loop over edges
    graph_edge = sites[i].edges;

  
    double v1[2] = {-1000,-1000};
    double v2[2] = {-1000,-1000};
    int count = 0;

    while(graph_edge){

      v2[0] = graph_edge->pos[1].x;
      v2[1] = graph_edge->pos[1].y;

      if ( is_equal_double(v1, v2, 2) != 1){
        voronoi[indx]->contour->vertex[count].x = v2[0];
        voronoi[indx]->contour->vertex[count].y = v2[1];
        ++count;

        v1[0] = v2[0];
        v1[1] = v2[1];
      }
      graph_edge = graph_edge->next;

    }


  }// end loop over sites


  // FILE *fp;
  // int numPoints = diagram.numsites;
  // char fileName[256];
  // for (int i = 0; i < numPoints; i++) {
  //    // write cell to file
  //     snprintf(fileName,sizeof(fileName),"Cells/%d.txt",i);
  //     // write clipped cell to files
  //     fp = fopen(fileName,"w");
  //     if ( fp != NULL)
  //       gpc_write_polygon(fp, 0, voronoi[i]);
  //     fclose(fp);
  
  // }

  jcv_diagram_free(&diagram);

  return voronoi;



  return 0;


}
