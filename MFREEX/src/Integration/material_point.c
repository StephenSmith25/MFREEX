#include "Integration/material_point.h"


#include "ShapeFunction/mls_shapefunction_materialpoint.h"
#include "Integration/DomainMaterialPoint.h"
#ifndef QUADRATURE_ORDER
#define QUADRATURE_ORDER 1
#endif


#define APPROX_NEIGHBOUR_NUMBER 15



static int IS_AXI = -1;



#ifndef DIM 
#define DIM 2
#endif 

// distance function: returns normalised distance from x to x_I
inline double distance_function(double *x_I, double * x, MAT * invMI)
{

	// norm_distance = beta * (x-x_I)' * invM_I * (x-x_I)
	double norm_distance = 0;

#if DIM == 2
	double x_m_xI[2] = {x[0] - x_I[0], 
						x[1] - x_I[1]}; 
	norm_distance =(
	  (x_m_xI[0])*(invMI->me[0][0]* x_m_xI[0] + invMI->me[0][1]* x_m_xI[1])	 // row 1 
	+ (x_m_xI[1])*(invMI->me[1][0]* x_m_xI[0] + invMI->me[1][1]* x_m_xI[1])  // row 2
	);

#elif DIM == 3
	double x_m_xI[3] = {x[0] - x_I[0],
						x[1] - x_I[1],
						x[2] - x_I[2]}; 
	norm_distance = (
	      (x_m_xI[0])*(invMI->me[0][0]* x_m_xI[0] + invMI->me[0][1]* x_m_xI[1] + invMI->me[0][2]* x_m_xI[2]) // row 1 
		+ (x_m_xI[1])*(invMI->me[1][0]* x_m_xI[0] + invMI->me[1][1]* x_m_xI[1] + invMI->me[1][2]* x_m_xI[2] ); // row 2
		+ (x_m_xI[1])*(invMI->me[2][0]* x_m_xI[0] + invMI->me[2][1]* x_m_xI[1] + invMI->me[2][2]* x_m_xI[2] ) // row 3
		); 
#endif

	return sqrt(norm_distance);

}

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


MATERIAL_POINTS * create_material_points(void * cells, 
	int IS_AXI_I, 
	int dim, 
	char * integration_type,
	MATERIAL_TYPE material_type, 
	CELLS * grid, 
	double rho,  
	double beta, 
	meshfreeDomain * mfree)

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


			MPS[i]->temperature = 0; 


			// Density
			MPS[i]->rho = rho;


			// Stress Voigt
			MPS[i]->stressVoigt = v_get(dim_p);

			if ( IS_AXI == 1)
			{
			double r = _scni_obj->scni[i]->r;
				MPS[i]->INTEGRATION_FACTOR = r*2*PI;
			}else{
				MPS[i]->INTEGRATION_FACTOR = 1;
			}

			MPS[i]->stateOld= new_material_state(MPS[i]->temperature, material_type,
				dim, IS_AXI);
			MPS[i]->stateNew = new_material_state(MPS[i]->temperature, material_type,
				dim, IS_AXI);

			MPS[i]->Jn = 1.00;



		}



		material_points->MP = MPS;
		material_points->dim = dim;
		material_points->IS_AXI = IS_AXI;




	}else if ( strcmp(integration_type,"TRIANGLE") == 0)

	{
		TRIANGLE * tri = (TRIANGLE *) cells;

		double * tri_points = tri->points;
		int * triangles = tri->triangles;
		int number_of_triangles = tri->num_triangles;

		
		int * quad_orders = malloc(number_of_triangles*sizeof(int));


		// quadrature order
		for ( int i = 0 ; i < number_of_triangles ; i++)
		{
			quad_orders[i] = QUADRATURE_ORDER;
			// if ( ( i == 21) || ( i == 77) || (i == 80) )
			// {
			// 	quad_orders[i] = 2;
			// }
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
			MPS[i]->neighbours = iv_get(APPROX_NEIGHBOUR_NUMBER);


			// Support parameters 
			MPS[i]->MI = m_get(dim,dim);
			MPS[i]->invMI = m_get(dim,dim);

			// set the initial domain of the material points 
			setDomainMaterialPoint(mfree->nodes, grid, MPS[i]);


			MPS[i]->shape_function = NULL;
			MPS[i]->shape_function = mls_shapefunction_materialpoint(MPS[i],2,mfree->nodes);



			// getting Bmat 
			// Strain Displacement relationship 
			MPS[i]->B = BMAT(MNULL,MPS[i]->shape_function,dim,IS_AXI,MPS[i]->coords_n_1[0]);


			// Matricies used in internal force calculation
			MPS[i]->fInt = v_get(dim*MPS[i]->num_neighbours);

			// Deformation gradient 
			MPS[i]->F_n = m_get(dim_s,dim_s);
			MPS[i]->inc_F = m_get(dim_s,dim_s);
			m_ident(MPS[i]->inc_F);
			m_ident(MPS[i]->F_n);
			MPS[i]->Jn = 1.00;

			// Temperature
			MPS[i]->temperature = 0;
			if ( mfree->temperatures != NULL){
				for ( int k = 0 ; k < MPS[i]->num_neighbours ; k++)
				{
					int index = MPS[i]->neighbours->ive[k];
					MPS[i]->temperature += MPS[i]->shape_function->phi->ve[k]*mfree->temperatures[index];
				}
			}


			// Density
			MPS[i]->rho = rho;


			MPS[i]->index = i;

			// Stress Voigt
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


			// Create material state for eahc material point 
			MPS[i]->stateOld= new_material_state(MPS[i]->temperature, material_type,
				dim, IS_AXI);
			MPS[i]->stateNew = new_material_state(MPS[i]->temperature, material_type,
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
//static call_count = 0;

int write_domains(char * filename, MATERIAL_POINTS * MPS)
{
	FILE * fp = fopen(filename,"w");


	int k = 0;
	for ( int k = 0 ; k < MPS->num_material_points ; k++ )
	{
		double M11 =  MPS->MP[k]->invMI->me[0][0];
		double M12 =  MPS->MP[k]->invMI->me[0][1];
		double M21 =  MPS->MP[k]->invMI->me[1][0];
		double M22 =  MPS->MP[k]->invMI->me[1][1]; 
		double theta = MPS->MP[k]->theta;
		// print domains
		fprintf(fp,"%lf,%lf,%lf,%lf,%lf\n",M11,M12,M21,M22,theta);


	}

	fclose(fp);

	// 	++call_count;


	// if ( call_count ==50)
	// {
	// 	exit(0);
	// }

	return 0;

}
int write_material_points(char * filename, MATERIAL_POINTS * MPS)
{
	FILE * fp = fopen(filename,"w");
	fprintf(fp,"x,y,z,s11,s22,s12\n");


	for ( int i = 0 ; i < MPS->num_material_points ; i++ )
	{
		double s11 = MPS->MP[i]->stateNew->sigma->me[0][0];
		double s22 = MPS->MP[i]->stateNew->sigma->me[1][1];
		double s12 = MPS->MP[i]->stateNew->sigma->me[0][1];

		// print stress as well
		fprintf(fp,"%lf,%lf,%lf,%lf,%lf,%lf\n",MPS->MP[i]->coords_n_1[0],MPS->MP[i]->coords_n_1[1],MPS->MP[i]->coords_n_1[2],
			s11,s22,s12
		);
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
	if ( IS_AXI == 1)
	{
		MP->inc_F->me[2][2] = 1.00;
	}
	double inc_Jacobian = determinant(MP->inc_F);

	// Save old volume 
	double int_factor_old = MP->INTEGRATION_FACTOR;
	double volume_old = MP->volume*int_factor_old;


	// update volume 
	MP->volume = MP->volume*inc_Jacobian;


	/* UPDATE MATERIAL POINT COORDINATES*/
	for ( int i = 0 ; i < MP->num_neighbours ; i++)
	{	
		for ( int k = 0 ; k < dim ; k++)
		{
			MP->coords_n[k] = MP->coords_n_1[k];
		}
	}

	if ( IS_AXI == 1)
	{
		MP->INTEGRATION_FACTOR = MP->coords_n_1[0]*2*PI;
	}else{
		// leave integration factor as 1
	}

	double volume_new = MP->volume*MP->INTEGRATION_FACTOR;
	MP->rho = MP->rho * (volume_old/volume_new);







	// update domain of material point 
	updateDomainMaterialPoint(NODES, grid, MP);

	// if ( MP->index == 1412)
	// {
	// 	//iv_foutput(stdout, MP->neighbours);
	// }


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





	return MP;
}

int RangeSearchMaterialPoint(
MATERIAL_POINT * MP, MAT * nodes,CELLS * cells)
{

	// NOTE - IMPLEMENTED JUST FOR 2D AT THE MOMENT d
	// Find all nodes that are within this range 
	int num_neighbours = 0;
	int dim = nodes->n;
	iv_zero(MP->neighbours);

	double * x = MP->coords_n_1;
	double * x_I;
	
	// convert x into its i,j components
	int j = (int)floor((x[0]-cells->Bx_min)/cells->CX);
	int i = (int)floor((x[1]-cells->By_min)/cells->CY);

	// number of cells in each direction
	int nx = cells->nx;
	int ny = cells->ny;



	NODELIST * current_p = NULL;
	double * node_check = NULL;
	double distance = 0;
	for ( int indx_j = j-1 ; indx_j <= j+1 ; indx_j++)
		for ( int indx_i = i-1 ; indx_i <= i+1 ; indx_i++)
	{


		if (( indx_j >= 0) && (indx_j < nx) && ( indx_i >= 0) && (indx_i < ny))
		{
			current_p = cells->cells[indx_i][indx_j].nodes;
			while ( current_p != NULL)
			{

				x_I = nodes->me[current_p->node_number];

				double norm_distance = distance_function(x_I, x, MP->invMI);

				if ( norm_distance <= 1.00)
				{
					MP->neighbours->ive[num_neighbours] = current_p->node_number;
					++num_neighbours;
				}


				current_p = current_p->next;

			}


		}else{
			// skip
		}


	}
	MP->num_neighbours = num_neighbours;

	return 0;
}