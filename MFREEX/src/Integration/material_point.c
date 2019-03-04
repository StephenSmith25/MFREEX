#include "Integration/material_point.h"


#define QUADRATURE_ORDER 2



static inline MATERIAL_POINT * create_material_point(double coords[3], double volume)
{

	MATERIAL_POINT * mp = malloc(1*sizeof(MATERIAL_POINT));


	mp->volume = volume;
	mp->coords = malloc(3*sizeof(double));

	mp->coords[0] = coords[0];
	mp->coords[1] = coords[1];
	mp->coords[2] = coords[2];


	return mp; 

}


MATERIAL_POINTS * create_material_points(void * cells, int IS_AXI, int dim, 
	char * integration_type, char * material, double rho,  meshfreeDomain * mfree)

{


	MATERIAL_POINTS * material_points = malloc(1*sizeof(MATERIAL_POINTS));



	if ( strcmp(integration_type,"SCNI") == 0)
	{
		voronoi_diagram * vor = (voronoi_diagram*) cells;

		int is_stabalised = 0;
		SCNI_OBJ * _scni_obj = NULL;

		_scni_obj = generate_scni(vor, NULL , is_stabalised, IS_AXI, dim, mfree);

		// assign memory for each material point
		MATERIAL_POINT ** MPS = malloc(vor->num_cells * sizeof(MATERIAL_POINT));


		material_points->sf_material_points = mls_shapefunction(mfree->nodes, 1, mfree);


		material_points->num_material_points = vor->num_cells;

		// create a material state for each point ( dependent on the material )
		double point_coords[3] = {0,0,0};


		for ( int i = 0 ; i < material_points->num_material_points ; i++)
		{

			// Generate the coordinates
			for  (  int k = 0 ; k < dim ; k++)
			{
				point_coords[k] = mfree->nodes->me[i][k];
			}

			// Create the material points
			MPS[i] = create_material_point(point_coords,_scni_obj->scni[i]->area);

			// Shape function matricies 
			MPS[i]->B = _scni_obj->scni[i]->B;
			MPS[i]->neighbours = _scni_obj->scni[i]->sfIndex;
			MPS[i]->fInt = _scni_obj->scni[i]->fInt;
			MPS[i]->F_r = _scni_obj->scni[i]->F_r;
			MPS[i]->rho = rho;

		}



		material_points->MP = MPS;
		material_points->dim = dim;
		material_points->IS_AXI = IS_AXI;
		material_points->MP_COORDS = mfree->nodes;




	}else if ( strcmp(integration_type,"TRIANGLE") == 0)

	{
		struct triangulateio * tri = (struct triangulateio *) cells;

		double * tri_points = tri->pointlist;
		int * triangles = tri->trianglelist;
		int number_of_triangles = tri->numberoftriangles;
		int * quad_orders = malloc(number_of_triangles*sizeof(int));


		// quadrature order
		for ( int i = 0 ; i < number_of_triangles ; i++)
		{
			quad_orders[i] = QUADRATURE_ORDER;
		}	

		//  create quadrature points
		QUAD_TRIANGLE * quad_triangles = create_triangle_quadrature_points(tri_points, 
		triangles,number_of_triangles, quad_orders );



		int number_of_quadrature_points = quad_triangles->QUAD_POINTS->m;
		material_points->num_material_points = number_of_quadrature_points;


		MATERIAL_POINT ** MPS = malloc(number_of_quadrature_points * sizeof(MATERIAL_POINT));



		double point_coords[3] = {0,0,0};


		// Create BMAT;
		shape_function_container * sf_material_points = mls_shapefunction(quad_triangles->QUAD_POINTS, 2, mfree);
		material_points->sf_material_points = sf_material_points;



		for ( int i = 0 ; i < number_of_quadrature_points ; i++)
		{	


			// Generate the coordinates
			for  (  int k = 0 ; k < dim ; k++)
			{
				point_coords[k] = quad_triangles->QUAD_POINTS->me[i][k];
			}

			// Create the material points
			MPS[i] = create_material_point(point_coords,quad_triangles->VOLUMES->ve[i]);

			// generate BMAT and neighbour list;
			MPS[i]->B = generate_Bmat(sf_material_points->sf_list[i]->dphi,dim,IS_AXI,-1);
			MPS[i]->neighbours = sf_material_points->sf_list[i]->neighbours;
			MPS[i]->phi = sf_material_points->sf_list[i]->phi;

			// Matricies used in internal force calculation
			MPS[i]->fInt = v_get(dim*MPS[i]->neighbours->max_dim);
			MPS[i]->F_r = m_get(dim,dim);
			m_ident(MPS[i]->F_r);
			MPS[i]->rho = rho;




		}




		material_points->MP_COORDS = quad_triangles->QUAD_POINTS;
		material_points->MP = MPS;
		material_points->dim = dim;
		material_points->IS_AXI = IS_AXI;





	}


	return material_points;
}

static inline int delete_material_point()
{
	return 0;
}


int write_material_points(char * filename, MATERIAL_POINTS * MPS)
{
	FILE * fp = fopen(filename,"w");


	for ( int i = 0 ; i < MPS->num_material_points ; i++ )
	{
		fprintf(fp,"%lf,%lf,%lf\n",MPS->MP[i]->coords[0],MPS->MP[i]->coords[1],MPS->MP[i]->coords[2]);
	}

	fclose(fp);


	return 0;
}


int set_material_point_temperature(int i , double temp);


int update_material_point();