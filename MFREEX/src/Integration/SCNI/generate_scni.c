#include "Integration/SCNI/generate_scni.h"
int isCCW(double * a, double * b, double * c){

	double det = (b[0]*c[1] + a[0]*b[1] + a[1]*c[0]) - (a[1]*b[0] + b[1]*c[0] + a[0]*c[1]);
	
	if ( det > 0 ){
		return 1;
	}else{
		return 0;
	}
}



// find if int a is in the list of integers stored in b
int findInt( int a, int * b, int size_b)
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


SCNI_OBJ * generate_scni(voronoi_diagram * voronoi, char * type, int is_stabalised, int is_AXI, int dim, meshfreeDomain * Mfree){

	// cell points and meshfree nodes
	MAT * cell_verticies = voronoi->verticies;
	MAT * nodes = Mfree->nodes;

	printf("getting shape function value at cell veritices \n ");
	shape_function_container * sf_verticies = mls_shapefunction(cell_verticies, "linear", "cubic", 2, 1, Mfree);




	// if axi will have to find shape function at nodes to add the unsmoothed component
	// at some nodes will have to find, compute could be given as a vector 
	if ( is_AXI  == 1)
	{
		shape_function_container * sf_nodes = mls_shapefunction(nodes, "linear" , "quartic", 2, 3, Mfree);

	}


	int num_cells = voronoi->num_cells;

	printf("Looping over cells and construction B \n ");
	
	VEC * phi = VNULL;
	IVEC * neighbours = IVNULL;


	// cell combined neighbours
	// appr
	int approx_cell_sf_index_length = 20;


	SCNI ** scni_  = malloc(num_cells*sizeof(SCNI*));

	int i;
	IVEC * temp_vec = iv_get(approx_cell_sf_index_length);

	for ( int k = 0 ; k < approx_cell_sf_index_length ; k++)
	{
		temp_vec->ive[k] = 1;
		// cell variables
	}
	int length_cell_index = 0;

	int * cell_index = NULL;
	int num_cell_verticies = 0;

		// segment normals and segment lengths
	double n_m[2];
	double n_m_1[2];
	double seg_length_m =0;
	double seg_length_m_1 = 0;

		// cell area and center
	double center[2] = {0,0};
	double area = 0;
		// size of sfIndex
	MAT * bI = m_get(approx_cell_sf_index_length,2);

	for ( i = 0 ; i < num_cells; i++)
	{
		scni_[i] = malloc(1*sizeof(SCNI));
		cell_index = voronoi->index[i];
		num_cell_verticies = voronoi->num_cell_verticies[i];

		// smoothed shape function gradients
		IVEC * cell_sf_index = iv_get(approx_cell_sf_index_length);
		iv_sub(cell_sf_index,temp_vec, cell_sf_index);
		bI = m_resize(bI,approx_cell_sf_index_length,2);
		m_zero(bI);



		// re initialise parameters
		area = 0;
		center[0] = 0;
		center[1] = 0;
		length_cell_index = 0;


		// check orientation of the polygon
		double * v1 = cell_verticies->me[cell_index[0]];
		double * v2 = cell_verticies->me[cell_index[1]];
		double * v3 = cell_verticies->me[cell_index[2]];

		int isCCW_ = isCCW(v1,v2,v3);
		int normalFactor;
		// needed to make sure the normal is outward pointing when doing contour integration
		if ( isCCW_ == 1){
			normalFactor = 1;
		}else{
			normalFactor = -1;
		}

		//
		seg_length_m = sqrt(pow(v2[0]-v1[0],2) + pow(v2[1]-v1[1],2));  
		n_m[0] = (normalFactor)*1.00*(v2[1] - v1[1])/seg_length_m;
		n_m[1] = (normalFactor)*-1.00*(v2[0] - v1[0])/seg_length_m;

		double n[num_cell_verticies][2];
		double l[num_cell_verticies];
		for ( int k = 0 ; k < num_cell_verticies ; k++)
		{
			v1 = cell_verticies->me[cell_index[k]];
			if ( k == num_cell_verticies -1)
			{
				v2 = cell_verticies->me[cell_index[0]];

			}else{
				v2 = cell_verticies->me[cell_index[k+1]];

			}

			// segment length
			l[k] = sqrt(pow(v2[0]-v1[0],2) + pow(v2[1]-v1[1],2));  
			// segment normals
			n[k][0] = (normalFactor)*1.00*(v2[1] - v1[1])/seg_length_m;
			n[k][1] = (normalFactor)*-1.00*(v2[0] - v1[0])/seg_length_m;



			area += (v1[0]*v2[1] - v2[0]*v1[1]);
			center[0] += v2[0];
			center[1] += v2[1];
		}

		center[0] = center[0]/num_cell_verticies;
		center[1] = center[1]/num_cell_verticies;
		area = 0.5*fabs(area);




		// construct smoothed shape function gradients in the cell
		for ( int k = 0 ; k < num_cell_verticies ; k++)
		{
			// node 
			int indx;


			if ( k == num_cell_verticies -1)
			{
				indx = cell_index[0];

			}else{
				indx = cell_index[k+1];
			}



			// get phi at this point 
			phi = sf_verticies->sf_list[indx]->phi;
			neighbours = sf_verticies->sf_list[indx]->neighbours;

			// recursive defitnion b_Ii = sum_Ns [ 0.5 * ( n_im * l_m + n_i(m+1) * l_(m+1) ) phi_I ( M+1)]
			for ( int j = 0 ; j < neighbours->max_dim ; j++)
			{
				// find if node_nr is already in 
				int node_nr = neighbours->ive[j];
				int idx = findInt(node_nr, cell_sf_index->ive, length_cell_index);

					// point was already in cell_sf_index
				if ( idx != -1)
				{
					// add b_Ii contrbution to b_I;
					bI->me[idx][0] += 0.5*(n[k][0]*l[k] + n[k+1][0]*l[k+1])*phi->ve[j];
					bI->me[idx][1] += 0.5*(n[k][1]*l[k] + n[k+1][1]*l[k+1])*phi->ve[j];

				}else{

					// increment counter
					++length_cell_index;

					// check sufficent memory allocated to write into cell_sf_index and bI
					if ( length_cell_index > cell_sf_index->max_dim )
					{	
						MAT * temp_m = m_copy(bI,MNULL);
						IVEC * temp_iv = iv_copy(cell_sf_index,IVNULL);
						m_resize(bI,length_cell_index+10,2);
						iv_resize(cell_sf_index, length_cell_index+10);
						m_move(temp_m, 0, 0,temp_m->m, temp_m->n, bI, 0, 0);
						iv_move(temp_iv, 0, temp_iv->max_dim, cell_sf_index, 0);
						M_FREE(temp_m);
						IV_FREE(temp_iv);

					}

					// add contribution to bI and cell_sf_index 
					cell_sf_index->ive[length_cell_index-1] = node_nr;
					bI->me[length_cell_index-1][0] += 0.5*(n[k][0]*l[k] + n[k+1][0]*l[k+1])*phi->ve[j];
					bI->me[length_cell_index-1][1] += 0.5*(n[k][1]*l[k] + n[k+1][1]*l[k+1])*phi->ve[j];

				}
			}
			// int factor = 0.5 * ( n_im * l_m + n_i(m+1) * l_(m+1) ) 



		} // end of loop over cell verticies


		// m_foutput(stdout, bI);
		// printf("got to resizing \n");
		// printf("length_cell_index = %d\n",length_cell_index);
		// printf("length of cell_sf_index = %d \n", cell_sf_index->max_dim );
		// printf("i = %d \n ", i);
		if ( cell_sf_index->max_dim > length_cell_index   )
		{	
			MAT * temp_m = m_copy(bI,MNULL);
			IVEC * temp_iv = iv_copy(cell_sf_index,IVNULL);
			m_resize(bI,length_cell_index,2);
			iv_resize(cell_sf_index, length_cell_index);
			m_move(temp_m, 0, 0,length_cell_index-1, temp_m->n, bI, 0, 0);
			iv_move(temp_iv, 0, length_cell_index-1, cell_sf_index, 0);
			M_FREE(temp_m);
			IV_FREE(temp_iv);

		}

		// generate Bmat and set up SCNI structure
		sm_mlt(1.000/area, bI, bI);
		scni_[i]->B = generate_Bmat(bI,dim,is_AXI,-1);
		scni_[i]->sfIndex = cell_sf_index;
		scni_[i]->area = area;
		scni_[i]->center = center;
		scni_[i]->fInt = v_get(dim*cell_sf_index->max_dim);

	} // end of loop over cells


	M_FREE(bI);


	printf("Finished looping over cells and constructing B \n ");


	free_shapefunction_container(sf_verticies);


	SCNI_OBJ * scni_obj = malloc(1*sizeof(SCNI_OBJ));
	scni_obj->scni = scni_;
	scni_obj->num_points = num_cells;

	return scni_obj;
}


int free_scni(SCNI * _scni){


	return 0;
}
