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
	shape_function_container * sf_nodes ;


	int dim_B = dim;

	// if axi will have to find shape function at nodes to add the unsmoothed component
	// at some nodes will have to find, compute could be given as a vector 
	if ( is_AXI  == 1)
	{
		sf_nodes = mls_shapefunction(nodes, "linear" , "cubic", 2, 2, Mfree);
		dim_B = 3;

	}

	int num_cells = voronoi->num_cells;
	
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

	// cell area and center
	double center[2] = {0,0};
	double area = 0;
		// size of sfIndex
	MAT * bI = m_get(Mfree->num_nodes,dim_B);
	double area_sum = 0;

	for ( i = 0 ; i < num_cells; i++)
	{
		scni_[i] = malloc(1*sizeof(SCNI));
		cell_index = voronoi->index[i];
		num_cell_verticies = voronoi->num_cell_verticies[i];

		// smoothed shape function gradients
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



		double n[num_cell_verticies][2];
		double l[num_cell_verticies];
		for ( int k = 0 ; k < num_cell_verticies ; k++)
		{
			v1 = cell_verticies->me[cell_index[k]];
			if ( k == (num_cell_verticies -1))
			{
				v2 = cell_verticies->me[cell_index[0]];

			}else{
				v2 = cell_verticies->me[cell_index[k+1]];

			}

			// segment length
			l[k] = sqrt(pow(v2[0]-v1[0],2) + pow(v2[1]-v1[1],2));  
			// segment normals
			n[k][0] = (normalFactor)*1.00*(v2[1] - v1[1])/l[k];
			n[k][1] = (normalFactor)*-1.00*(v2[0] - v1[0])/l[k];
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
			int m_1 = 0;
			int m_2 = 0;


			if ( k == num_cell_verticies -1)
			{
				indx = cell_index[0];
				m_1 = k;
				m_2 = 0;

			}else{
				indx = cell_index[k+1];
				m_1 = k;
				m_2 = k+1;
			}



			// get phi at this point 
			phi = sf_verticies->sf_list[indx]->phi;
			neighbours = sf_verticies->sf_list[indx]->neighbours;


			// recursive defitnion b_Ii = sum_Ns [ 0.5 * ( n_im * l_m + n_i(m+1) * l_(m+1) ) phi_I ( M+1)]
			for ( int j = 0 ; j < neighbours->max_dim ; j++)
			{
				int indx = neighbours->ive[j];
				// add contribution to bI and cell_sf_index 
				bI->me[indx][0] += 0.5*(n[m_1][0]*l[m_1] + n[m_2][0]*l[m_2])*phi->ve[j];
				bI->me[indx][1] += 0.5*(n[m_1][1]*l[m_1] + n[m_2][1]*l[m_2])*phi->ve[j];

			}


		} // end of loop over cell verticies
		int length_index = 0;

		for ( int k = 0 ; k < Mfree->num_nodes; k++)
		{
			if (( bI->me[k][0] != 0) || ( bI->me[k][1] !=0 ))
			{
				++length_index;
			}

		}

		IVEC * cell_sf_index = iv_get(length_index); 
		MAT * bI_n = m_get(length_index,dim_B);
		length_index = 0;
		for ( int k = 0 ; k < Mfree->num_nodes; k++)
		{
			if (( bI->me[k][0] != 0) || ( bI->me[k][1] !=0 ))
			{
				cell_sf_index->ive[length_index] = k;
				bI_n->me[length_index] = bI->me[k];  
				++length_index;
			}

		}
		sm_mlt(1.000/area, bI_n, bI_n);

		area_sum += area;

		if ( is_AXI == 1)
		{
			neighbours = sf_nodes->sf_list[i]->neighbours;
			double r = nodes->me[i][0];

			for ( int k = 0 ; k < neighbours->max_dim ; k++)
			{
				//double r = nodes->me[neighbours->ive[k]][0];
				if ( r > 0)
				{
					int indx = findInt(neighbours->ive[k], cell_sf_index->ive, cell_sf_index->max_dim);
					bI_n->me[indx][2] += sf_nodes->sf_list[i]->phi->ve[k]/r;
				}
				else{
					int indx = findInt(neighbours->ive[k], cell_sf_index->ive, cell_sf_index->max_dim);
					bI_n->me[indx][2] += sf_nodes->sf_list[i]->dphi->me[k][0];
				}


			}
			
			
		}


		

		// generate Bmat and set up SCNI structure
		scni_[i]->B = generate_Bmat(bI_n,dim,is_AXI,-1);
		M_FREE(bI_n);
		scni_[i]->sfIndex = cell_sf_index;
		scni_[i]->area = area;
		scni_[i]->r = center[0];

		scni_[i]->fInt = v_get(dim*cell_sf_index->max_dim);
		
		if ( is_AXI == 1){
			scni_[i]->F_r = m_get(3,3);
		}else{
			scni_[i]->F_r = m_get(dim,dim);
		}
		m_ident(scni_[i]->F_r);

	} // end of loop over cells


	M_FREE(bI);




	free_shapefunction_container(sf_verticies);

	if ( is_AXI == 1)
	{
		free_shapefunction_container(sf_nodes);

	}


	SCNI_OBJ * scni_obj = malloc(1*sizeof(SCNI_OBJ));
	scni_obj->scni = scni_;
	scni_obj->num_points = num_cells;
	printf("area = %lf \n",area_sum);
	return scni_obj;
}


int free_scni(SCNI * _scni){


	return 0;
}
