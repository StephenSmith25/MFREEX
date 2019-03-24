#include "Integration/material_point.h"


#include "ShapeFunction/mls_shapefunction_materialpoint.h"
#include "Integration/DomainMaterialPoint.h"
#ifndef QUADRATURE_ORDER
#define QUADRATURE_ORDER 1
#endif



static int IS_AXI = -1;



static inline MATERIAL_POINT * create_material_point(double coords[3], double volume)
{

	MATERIAL_POINT * mp = malloc(1*sizeof(MATERIAL_POINT));


	mp->volume = volume;
	mp->coords_n_1 = malloc(3*sizeof(double));
	mp->coords_n = malloc(3*sizeof(double));


	// Generate the coordinates
	for  (  int k = 0 ; k < 3 ; k++)
	{
		mp->coords_n_1[k] = coords[k];
		mp->coords_n[k] = coords[k];
	}



	return mp; 

}


MATERIAL_POINTS * create_material_points(void * cells, int IS_AXI_I, int dim, 
	char * integration_type, MATERIAL_TYPE material_type, 
	CELLS * grid, double rho,  double beta, meshfreeDomain * mfree)

{

	MATERIAL_POINTS * material_points = malloc(1*sizeof(MATERIAL_POINTS));

	int dim_s = dim;
	int dim_p = dim*dim;

	IS_AXI = IS_AXI_I;

	if ( IS_AXI == 1)
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


		material_points->num_material_points = vor->num_cells;

		// create a material state for each point ( dependent on the material )
		double point_coords[3] = {0,0,0};


		for ( int i = 0 ; i < material_points->num_material_points ; i++)
		{


			MPS[i]->coords_n_1 = malloc(dim*sizeof(double));
			MPS[i]->coords_n = malloc(dim*sizeof(double));


			// Generate the coordinates
			for  (  int k = 0 ; k < dim ; k++)
			{
				MPS[i]->coords_n_1[k] = mfree->nodes->me[i][k];
				MPS[i]->coords_n[k] = mfree->nodes->me[i][k];
			}

			// Create the material points
			MPS[i] = create_material_point(point_coords,_scni_obj->scni[i]->area);

			// Shape function matricies 
			MPS[i]->B = _scni_obj->scni[i]->B;
			//MPS[i]->shape_function = MLS_SHAPEFUNCTION(NULL,point_coords,1,mfree);
			MPS[i]->fInt = _scni_obj->scni[i]->fInt;

			// deformation gradients
			MPS[i]->inc_F = m_get(dim_s,dim_s);
			MPS[i]->F_n = _scni_obj->scni[i]->F_r;

			m_ident(MPS[i]->inc_F);
			m_ident(MPS[i]->F_n);




			MPS[i]->rho = rho;

			MPS[i]->stressVoigt = v_get(dim_p);

			if ( IS_AXI == 1)
			{
			double r = _scni_obj->scni[i]->r;
				MPS[i]->INTEGRATION_FACTOR = r*2*PI;
			}else{
				MPS[i]->INTEGRATION_FACTOR = 1;
			}


			MPS[i]->stateOld= new_material_state(0.00, material_type,
				dim, IS_AXI);
			MPS[i]->stateNew = new_material_state(0.00, material_type,
				dim, IS_AXI);

			MPS[i]->Jn = 1.00;



		}



		material_points->MP = MPS;
		material_points->dim = dim;
		material_points->IS_AXI = IS_AXI;




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




		for ( int i = 0 ; i < number_of_quadrature_points ; i++)
		{	

			// Generate the coordinates
			for  (  int k = 0 ; k < dim ; k++)
			{
				point_coords[k] = quad_triangles->QUAD_POINTS->me[i][k];
			}

			// Create the material points
			MPS[i] = create_material_point(point_coords,quad_triangles->VOLUMES->ve[i]);

			MPS[i]->beta = beta;
			MPS[i]->kernel_support=RADIAL;
			MPS[i]->neighbours = iv_get(50);

			// set the initial domain of the material points 
			setDomainMaterialPoint(mfree->nodes, grid, MPS[i]);



			MPS[i]->shape_function = NULL;
			MPS[i]->shape_function = mls_shapefunction_materialpoint(MPS[i],2,mfree->nodes);


			// generate BMAT 
			MPS[i]->B = BMAT(MNULL,MPS[i]->shape_function,dim,IS_AXI,MPS[i]->coords_n_1[0]);


			// Matricies used in internal force calculation
			MPS[i]->fInt = v_get(dim*MPS[i]->num_neighbours);

			// deformation gradients
			MPS[i]->F_n = m_get(dim_s,dim_s);
			MPS[i]->inc_F = m_get(dim_s,dim_s);
			m_ident(MPS[i]->inc_F);
			m_ident(MPS[i]->F_n);


			MPS[i]->Jn = 1.00;

			MPS[i]->rho = rho;
			MPS[i]->stressVoigt = v_get(dim_p);

			MPS[i]->temp = m_get(dim_s,dim_s);
			MPS[i]->temp_1 = m_get(dim_s,dim_s);

			// Modify integration factor if problem is axisymmetric
			if ( IS_AXI == 1)
			{
				MPS[i]->INTEGRATION_FACTOR = 2*PI*MPS[i]->coords_n_1[0];
			}else{
				MPS[i]->INTEGRATION_FACTOR = 1.00;
			}
			MPS[i]->stateOld= new_material_state(0.00, material_type,
				dim, IS_AXI);
			MPS[i]->stateNew = new_material_state(0.00, material_type,
				dim, IS_AXI);


		}



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
		fprintf(fp,"%lf,%lf,%lf\n",MPS->MP[i]->coords_n_1[0],MPS->MP[i]->coords_n_1[1],MPS->MP[i]->coords_n_1[2]);
	}

	fclose(fp);


	return 0;
}


void set_material_point_temperature(double temperature , MATERIAL_POINT * MP)
{
	MP->temperature = temperature;
}



MATERIAL_POINT * update_material_point(MATERIAL_POINT * MP, CELLS * grid, MAT * NODES, VEC * nodal_mass)
{


	// // update density material point
	int dim = NODES->n;
	int index = 0;



	// // Store F_n at each material point 
	MP->F_n = m_copy(MP->stateNew->F,MP->F_n);
	MP->Jn = determinant(MP->F_n);


	// Incremental deformation gradient
	double inc_Jacobian = determinant(MP->inc_F);

	// update volume and density of material points
	MP->rho = MP->rho / inc_Jacobian;
	MP->volume = MP->volume*inc_Jacobian;


	/* UPDATE MATERIAL POINT COORDINATES*/
	for ( int i = 0 ; i < MP->num_neighbours ; i++)
	{	
		index = MP->neighbours->ive[i];

		for ( int k = 0 ; k < dim ; k++)
		{
			if ( i == 0)
			{
				MP->coords_n[k] = MP->coords_n_1[k];
				MP->coords_n_1[k] = 0;	
			}

			MP->coords_n_1[k] += MP->shape_function->phi->ve[i] * NODES->me[index][k];
		}
	}


	// m_inverse_small(MP->inc_F, MP->temp);

	// mtrm_mlt(MP->temp, MP->invMI, MP->temp_1);

	// m_zero(MP->invMI);
	// m_mlt(MP->temp_1,MP->temp,MP->invMI);


	updateDomainMaterialPoint(NODES, grid, MP);



	//REFORM SHAPE FUNCTIONS
	MP->shape_function = mls_shapefunction_materialpoint(MP, 2 , NODES);

	//Reform B matrix 
	MP->B = BMAT(MP->B,MP->shape_function,dim,IS_AXI,MP->coords_n_1[0]);





	// update mass vector
	for ( int k = 0 ; k < MP->num_neighbours ; k++)
	{
		index = MP->neighbours->ive[k];		
		nodal_mass->ve[index] += MP->INTEGRATION_FACTOR*MP->volume * 
		MP->rho*MP->shape_function->phi->ve[k];
	}

	MP->fInt = v_resize(MP->fInt, dim*MP->num_neighbours);

	if ( IS_AXI == 1)
	{
		MP->INTEGRATION_FACTOR = MP->coords_n_1[0]*2*PI;
	}else{
		// leave integration factor as 1
	}




	return MP;
}