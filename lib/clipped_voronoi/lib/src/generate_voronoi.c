#include "generate_voronoi.h"
#include "clipVoronoi.h"
#include "callVoronoi.h"
#include <math.h>


// find if the point is already a member of the vertex list
static inline int is_member(double * point, double ** points, int numPoints, int dim ){


	double tol_equal = 1e-5;

	int is_member_row = 0;
	int i = 0;


	while ((is_member_row != 1) && (i < numPoints)){

		if ( fabs(point[0] - points[i][0]) < tol_equal){
			int k = 1;
			is_member_row = 1;
			while (( k < dim) && ( is_member_row != 0))
			{	
				if ( fabs(point[k] - points[i][k]) < tol_equal){
					// could be a membmer of row i
					is_member_row = 1;
				}else{
					// point and row i are not equal
					is_member_row = 0;
				}
				// increment k
				++k;

			}// end of while
		}// end of if points are equal


		// increment i
		++i;
	}
	if ( is_member_row == 1)
	{
		return i-1;

	}else{
		return -1;
	}
}

voronoi_diagram * generate_voronoi(double *points, int * boundary, int numPoints, int numBoundary, int dim){





	// generate clipped voronoi diagram
	gpc_polygon ** voronoi_1 = callVoronoi(points,numPoints);
	clipVoronoi(&voronoi_1,points,boundary,numPoints,numBoundary);




	// move polygon into new voronoi structure
	MAT * vertex_list = m_get(numPoints,dim);
	m_ones(vertex_list);
	sm_mlt(1000, vertex_list, vertex_list);

	int count_vertex_list = 0;


	// 
	int total_verticies = 0;


	// set up output verticies
	voronoi_diagram * vor_out = malloc(sizeof(voronoi_diagram));
	vor_out->index = (int**)malloc(numPoints * sizeof(int*));
	vor_out->num_cell_verticies = malloc(numPoints*sizeof(int));
	vor_out->num_cells = numPoints;


	gpc_polygon * vor;

	for ( int i = 0 ; i < numPoints ; i++)
	{	

		// loop over each vertex of the cell
		vor = voronoi_1[i];
		int num_vertices = vor->contour->num_vertices;

		total_verticies += num_vertices;

		vor_out->index[i] = (int*) malloc(num_vertices * sizeof(int));
		vor_out->num_cell_verticies[i] = num_vertices;
		for ( int k = 0 ; k < num_vertices ; ++k ){
			
			double point[2] = {vor->contour->vertex[k].x,vor->contour->vertex[k].y};
			int index = -1;
			index = is_member(point, vertex_list->me, count_vertex_list, dim);


			if ( index != -1 )
			{
				vor_out->index[i][k] = index;

			}else{
				// add point to vertex list 
				// check storage is big enough
				int size_m = vertex_list->m;
				if ( count_vertex_list >= size_m )
				{

					MAT * temp = m_copy(vertex_list,MNULL);
					vertex_list = m_resize(vertex_list,count_vertex_list+50, dim);
					m_ones(vertex_list);
					sm_mlt(1000, vertex_list, vertex_list);
					vertex_list = m_move(temp,0,0,temp->m,temp->n,vertex_list,0,0);

					M_FREE(temp);

					for ( int j = 0 ; j < dim ; j++)
					{
						vertex_list->me[count_vertex_list][j] = point[j];

					}
					vor_out->index[i][k] = count_vertex_list;


				}else{
					for ( int j = 0 ; j < dim ; j++)
					{
						vertex_list->me[count_vertex_list][j] = point[j];

					}
					vor_out->index[i][k] = count_vertex_list;

				}
				// increment vertex list point counter
				count_vertex_list = count_vertex_list + 1;

			}// end of else

		}// end of loop over verticies
		gpc_free_polygon(voronoi_1[i]);
		free(voronoi_1[i]);
	} // end of loop over points

	free(voronoi_1);

	if ( vertex_list->m > count_vertex_list)
	{
		MAT * temp = m_copy(vertex_list,MNULL);
		m_resize(vertex_list,count_vertex_list, dim);
		vertex_list = m_move(temp,0,0,count_vertex_list,temp->n,vertex_list,0,0);
		M_FREE(temp);
	}

	vor_out->verticies = vertex_list;
	return vor_out;
}

int print_voronoi_diagram(FILE *fp, voronoi_diagram * vor)
{
	fprintf(fp,"%d %d\n",vor->verticies->m,2);

	for ( int i = 0 ; i < vor->verticies->m ; i++)
	{
		fprintf(fp, "%lf %lf \n",vor->verticies->me[i][0], vor->verticies->me[i][1] );
	}

	fprintf(fp, "%d\n", vor->num_cells );

	for ( int i = 0 ; i < vor->num_cells ; i++)
	{
		fprintf(fp, "%d ",vor->num_cell_verticies[i]);

		for ( int k = 0 ; k < vor->num_cell_verticies[i] ; k++)
		{
			fprintf(fp, "%d ",vor->index[i][k]);
		}
		fprintf(fp,"\n");
	}


	return 0;
}
