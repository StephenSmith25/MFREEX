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
	char * integration_type, MATERIAL_TYPE material_type, double rho,  meshfreeDomain * mfree)

{


	MATERIAL_POINTS * material_points = malloc(1*sizeof(MATERIAL_POINTS));

	int dim_s = dim;
	int dim_p = dim*dim;

	if ( IS_AXI)
	{
		dim_s = dim+1;
		dim_p = dim_p + 1;
	}


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
			MPS[i]->F_n = _scni_obj->scni[i]->F_r;
			MPS[i]->rho = rho;

			MPS[i]->stressVoigt = v_get(dim_p);

			if ( IS_AXI == 1)
			{
			double r = _scni_obj->scni[i]->r;
			MPS[i]->volume = MPS[i]->volume*r*2*PI;
			}


			MPS[i]->stateOld= new_material_state(0.00, material_type,
				dim, IS_AXI);
			MPS[i]->stateNew = new_material_state(0.00, material_type,
				dim, IS_AXI);




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
			MPS[i]->F_n = m_get(dim_s,dim_s);
			m_ident(MPS[i]->F_n);
			MPS[i]->rho = rho;
			MPS[i]->stressVoigt = v_get(dim_p);

			if ( IS_AXI == 1)
			{
				MPS[i]->volume = MPS[i]->volume*2*PI*MPS[i]->coords[0];
			}
			MPS[i]->stateOld= new_material_state(0.00, material_type,
				dim, IS_AXI);
			MPS[i]->stateNew = new_material_state(0.00, material_type,
				dim, IS_AXI);


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



int update_material_point(MATERIAL_POINT * MP, MAT * NODES)
{


	// update density material point
	int dim = NODES->n;

	// Store F_n at each material point 
	MP->F_n = MP->stateNew->F;

	// Incremental deformation gradient
	double inc_Jacobian = determinant(MP->inc_F);




	// update volume and density of material points
	MP->rho = MP->rho / inc_Jacobian;
	MP->volume = MP->volume*inc_Jacobian;


	/* UPDATE MATERIAL POINT COORDINATES*/
	for ( int i = 0 ; i < MP->neighbours->max_dim ; i++)
	{	
		int index = MP->neighbours->ive[i];

		for ( int k = 0 ; k < dim ; k++)
		{
			MP->coords[k] = MP->phi->ve[i] * NODES->me[index][k];
		}
	}


	// REFORM SHAPE FUNCTIONS
	






	// reform mass matrix ? maybe




	// update material point coordinates




	// incremental deformation gradient

	// total deformation gradient

	// update material point density

	// update material point volume



	// update shape functions
		// Error : Insuffcient nodal connectivity
	// 

	// Update B







	return 0;
}