#include "Integration/SCNI/generate_mscni.h"


static inline int isCCW(double * a, double * b, double * c){

	double det = (b[0]*c[1] + a[0]*b[1] + a[1]*c[0]) - (a[1]*b[0] + b[1]*c[0] + a[0]*c[1]);
	
	if ( det > 0 ){
		return 1;
	}else{
		return 0;
	}
}



// find if int a is in the list of integers stored in b
static inline int findInt( int a, int * b, int size_b)
{
	int i = 0;
	int b_temp = -1;

	for ( int i = 0 ; i < size_b ; i++)
	{
		int b_temp = b[i];


		if ( b_temp == a)
		{
			return i;
		}
	}

	return -1;
}


static inline double tri_area(double *a, double *b, double * c)
{
	return fabs( 0.5*(   ( a[0]*(b[1]-c[1]) +  b[0]*(c[1]-a[1]) + c[0]*(a[1]-b[1])) ));
}



MSCNI_OBJ * generate_mscni(voronoi_diagram * voronoi, char * type, int is_stabalised, int is_AXI, meshfreeDomain * Mfree)
{



	// cell points and meshfree nodes
	MAT * cell_verticies = voronoi->verticies;
	MAT * nodes = Mfree->nodes;

	int dim  = Mfree->dim;



	// get shape function at verticies
	shape_function_container * sf_verticies = mls_shapefunction(cell_verticies, 1, Mfree);
	int num_cells = voronoi->num_cells;

	MSCNI ** mscni_  = malloc(num_cells*sizeof(MSCNI*));
	MAT * centers = m_get(num_cells,2);
	double sum_x =0;
	double sum_y = 0;



	for ( int k = 0 ; k < num_cells; k++)
	{	
		sum_x = 0;
		sum_y = 0;
		int * cell_index = voronoi->index[k];
		int num_cell_verticies = voronoi->num_cell_verticies[k];		

		for ( int m = 0 ; m < num_cell_verticies ; m++)
		{
			double * v1 = cell_verticies->me[cell_index[m]];
			sum_x += v1[0];
			sum_y += v1[1];

		}
		centers->me[k][0] = sum_x/num_cell_verticies;
		centers->me[k][1] = sum_y/num_cell_verticies;

	}


	shape_function_container * sf_nodes = mls_shapefunction(centers, 1, Mfree);
	omp_set_num_threads(4);

#pragma omp parallel
{



	VEC * phi = VNULL;
	IVEC * neighbours = IVNULL;

	// cell variables
	int length_cell_index = 0;
	int * cell_index = NULL;
	int num_cell_verticies = 0;


	// cell area and center
	double center[2] = {0,0};
	double area =0;
	double n[3][2];
	double l[3];
	int tri_index[3];

	// shape function gradient matrix
	MAT * bI = m_get(Mfree->num_nodes,dim*8);


	// pointers to verticies
	double * v1;
	double * v2;
	double * v3;

#pragma omp for
	// construct MSCNI for each cell. 
	for ( int i = 0 ; i < num_cells ; i++)
	{

		// pointer to the index of verticies that make up cell i
		cell_index = voronoi->index[i];
		mscni_[i] = malloc(1*sizeof(MSCNI));
		num_cell_verticies = voronoi->num_cell_verticies[i];
		mscni_[i]->stabalised_volumes = malloc(num_cell_verticies*sizeof(double));
		mscni_[i]->stabalised_B = malloc(num_cell_verticies*sizeof(MAT*));

		// smmothed shape function gradient
		m_zero(bI);

		if ( dim*num_cell_verticies > bI->n)
		{
			m_resize(bI,bI->m,dim*num_cell_verticies);
		}
		// re initialise paraemters
		center[0] = 0 ; center[1] = 0;
		length_cell_index = 0;
		area = 0;


		// check orientation of the polygon
		v1 = cell_verticies->me[cell_index[0]];
		v2 = cell_verticies->me[cell_index[1]];
		v3 = cell_verticies->me[cell_index[2]];

		int isCCW_ = isCCW(v1,v2,v3);
		int normalFactor;
		// needed to make sure the normal is outward pointing when doing contour integration
		if ( isCCW_ == 1){
			normalFactor = 1;
		}else{
			normalFactor = -1;
		}
		// Loop over each sub cell
		for ( int k = 0 ; k < num_cell_verticies ; k++){

			// verticies that make the triangle
			tri_index[0] = cell_index[k];
			tri_index[2] = i;

			if ( k == num_cell_verticies -1)
			{
				tri_index[1] = cell_index[0];

			}else{
				tri_index[1] = cell_index[k+1];
			}



			v1 = cell_verticies->me[tri_index[0]];
			v2 = cell_verticies->me[tri_index[1]];
			v3 = centers->me[i];



			// triangle segment normals
			l[0] = sqrt(pow(v2[0]-v1[0],2) + pow(v2[1]-v1[1],2)); // v2-v1
			l[1] = sqrt(pow(v3[0]-v2[0],2) + pow(v3[1]-v2[1],2)); // v3-v2
			l[2] = sqrt(pow(v1[0]-v3[0],2) + pow(v1[1]-v3[1],2)); // v1-v3
			n[0][0] = (normalFactor)*1.00*(v2[1] - v1[1])/l[0];
			n[0][1] = (normalFactor)*-1.00*(v2[0] - v1[0])/l[0];
			n[1][0] = (normalFactor)*1.00*(v3[1] - v2[1])/l[1];
			n[1][1] = (normalFactor)*-1.00*(v3[0] - v2[0])/l[1];
			n[2][0] = (normalFactor)*1.00*(v1[1] - v3[1])/l[2];
			n[2][1] = (normalFactor)*-1.00*(v1[0] - v3[0])/l[2];


			// area of triangle
			mscni_[i]->stabalised_volumes[k] = tri_area(v1,v2,v3);

	


			/*--------------------------------------------------------------------------------------------*/
			/*--------------------------------------------------------------------------------------------*/
			/* recursive defitnion b_Ii = sum_Ns [ 0.5 * ( n_im * l_m + n_i(m+1) * l_(m+1) ) phi_I ( M+1)]*/ 
			/*--------------------------------------------------------------------------------------------*/
			/*--------------------------------------------------------------------------------------------*/
			
			phi = sf_verticies->sf_list[tri_index[1]]->phi;
			neighbours = sf_verticies->sf_list[tri_index[1]]->neighbours;


			// First segment
			for ( int j = 0 ; j < neighbours->max_dim ; j++)
			{
				int indx = neighbours->ive[j];
				// add contribution to bI and cell_sf_index 
				bI->me[indx][2*k] += 0.5*(n[0][0]*l[0] + n[1][0]*l[1])*phi->ve[j];
				bI->me[indx][2*k+1] += 0.5*(n[0][1]*l[0] + n[1][1]*l[1])*phi->ve[j];

			}	

			// second segment
			phi = sf_nodes->sf_list[i]->phi;
			neighbours = sf_nodes->sf_list[i]->neighbours;

			for ( int j = 0 ; j < neighbours->max_dim ; j++)
			{
				int indx = neighbours->ive[j];
				// add contribution to bI and cell_sf_index 
				bI->me[indx][2*k] += 0.5*(n[1][0]*l[1] + n[2][0]*l[2])*phi->ve[j];
				bI->me[indx][2*k+1] += 0.5*(n[1][1]*l[1] + n[2][1]*l[2])*phi->ve[j];

			}	

			// third segment
			phi = sf_verticies->sf_list[tri_index[0]]->phi;
			neighbours = sf_verticies->sf_list[tri_index[0]]->neighbours;

			for ( int j = 0 ; j < neighbours->max_dim ; j++)
			{
				int indx = neighbours->ive[j];
				// add contribution to bI and cell_sf_index 
				bI->me[indx][2*k] += 0.5*(n[2][0]*l[2] + n[0][0]*l[0])*phi->ve[j];
				bI->me[indx][2*k+1] += 0.5*(n[2][1]*l[2] + n[0][1]*l[0])*phi->ve[j];

			}	





		


		} // end of loop over sub cells


		IVEC * cell_sf_index = iv_get(21);
		int count = 0;
		for ( int k = 0 ; k < Mfree->num_nodes; k++)
		{
			int j = 0;
			int test = 0;
			while ((j < dim*num_cell_verticies) && (test == 0)){
					// check if row entry is non zero
				if ( fabs(bI->me[k][j]) > 1e-14){
					test = 1;

					if ( count >= cell_sf_index->max_dim)
					{
						IVEC * temp = iv_copy(cell_sf_index,IVNULL);
						iv_resize(cell_sf_index, count+5);
						cell_sf_index->max_dim = count+5;
						iv_move(temp, 0, temp->max_dim, cell_sf_index, 0);
						IV_FREE(temp);
					}
					cell_sf_index->ive[count] = k;
					++count; 
				}
					// else move to next entry
				++j;
			}
		}

		// resize cell_sf_index
		if ( cell_sf_index->max_dim > count)
		{
			IVEC * temp = iv_copy(cell_sf_index,IVNULL);
			iv_resize(cell_sf_index, count);
			cell_sf_index->max_dim = count;
			iv_move(temp, 0, count, cell_sf_index, 0);
			IV_FREE(temp);
		}

		// shape function gradient matricies
		int length_index = cell_sf_index->max_dim;
		MAT * bI_n = m_get(length_index,dim*num_cell_verticies);



		for ( int k = 0 ; k < cell_sf_index->max_dim; k++)
		{
			int indx = cell_sf_index->ive[k];
			for ( int j = 0 ; j < num_cell_verticies ; j++)
			{
				double volume = mscni_[i]->stabalised_volumes[j];
				if ( volume > 0){
				if ( dim ==2)
				{
					bI_n->me[k][2*j] = (1.00/volume)* bI->me[indx][2*j];
					bI_n->me[k][2*j+1] = (1.00/volume)* bI->me[indx][2*j+1];


				}
				}

			}

		}


		// B matrix 
		double num_r;
		double num_c;
		if ( is_AXI == 1)
		{
			num_r = 5;
		}else{
			num_r = dim*dim;
		}
		num_c = dim*length_index;


		// generate Bmat and set up SCNI structure
		double volume = 0;
		mscni_[i]->B = m_get(num_r,num_c);
		MAT * B_temp = m_get(length_index,dim);

		for ( int k = 0 ; k < num_cell_verticies ; k++)
		{	
			m_move(bI_n, 0, 2*k, length_index, 2, B_temp, 0, 0);

			mscni_[i]->stabalised_B[k] = generate_Bmat(B_temp,dim,is_AXI,-1);


			//ms_mltadd(mscni_[i]->B, mscni_[i]->stabalised_B[k], mscni_[i]->stabalised_volumes[k],mscni_[i]->B);
			volume += mscni_[i]->stabalised_volumes[k];
			__mltadd__(mscni_[i]->B->base, mscni_[i]->stabalised_B[k]->base, mscni_[i]->stabalised_volumes[k],num_r*num_c);


		}


		__smlt__(mscni_[i]->B->base, 1.000/volume, mscni_[i]->B->base, num_r*num_c);

		for ( int k = 0 ; k < num_cell_verticies ; k++)
		{
			__sub__(mscni_[i]->B->base, mscni_[i]->stabalised_B[k]->base,mscni_[i]->stabalised_B[k]->base , num_r*num_c);
		}

		mscni_[i]->sfIndex = cell_sf_index;
		mscni_[i]->area = volume;
		mscni_[i]->center = NULL;
		mscni_[i]->fInt = v_get(dim*cell_sf_index->max_dim);
		mscni_[i]->num_sub_cells = num_cell_verticies;

		if ( is_AXI == 1){
			mscni_[i]->F_r = m_get(3,3);
		}else{
			mscni_[i]->F_r = m_get(dim,dim);
		}
		m_ident(mscni_[i]->F_r);



		M_FREE(B_temp);
		M_FREE(bI_n);

	} // end of loop over cells

	M_FREE(bI);
}
	free_shapefunction_container(sf_verticies);
	free_shapefunction_container(sf_nodes);
	M_FREE(centers);

	MSCNI_OBJ * scni_obj = malloc(1*sizeof(SCNI_OBJ));
	scni_obj->scni = mscni_;
	scni_obj->num_points = num_cells;

	return scni_obj;
}