#include "Integration/SCNI/generate_scni.h"

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


SCNI ** generate_scni(voronoi_diagram * voronoi, char * type, int is_stabalised, int is_AXI, int dim, meshfreeDomain * Mfree){

	/*  Point of shape function construction */
	MAT * point = m_get(1,2);
	/*  Shape function structure */


	MAT * cell_verticies = voronoi->verticies;
	MAT * nodes = Mfree->nodes;
	
	printf("getting shape function value at cell veritices \n ");
	shape_function_container * sf_verticies = mls_shapefunction(cell_verticies, "quadratic", "quartic", 2, 1, Mfree);


	// if axi will have to find shape function at nodes to add the unsmoothed component
	// at some nodes will have to find, compute could be given as a vector 
	if ( is_AXI  == 1)
	{
		shape_function_container * sf_nodes = mls_shapefunction(nodes, "quadratic" , "quartic", 2, 3, Mfree);
	
	}


	int num_cells = voronoi->num_cells;

	printf("Looping over cells and construction B \n ");
	
	VEC * phi = VNULL;
	IVEC * neighbours = IVNULL;


	// cell combined neighbours
	// appr
	int approx_cell_sf_index_length = 20;
	IVEC * cell_sf_index = iv_get(approx_cell_sf_index_length);
	// smoothed shape function gradients
	MAT * bI = m_get(approx_cell_sf_index_length,2);


	// cell variables
	int length_cell_index = 0;
	int * cell_index = NULL;
	int num_cell_verticies = 0;

	// segment normals and segment lengths
	double n_m[2];
	double n_m_1[2];
	double seg_length_m =0;
	double seg_length_m_1 = 0;

	// cell area and center
	double area = 0;
	double center[2] = {0,0};
	double total_area = 0;



	// loop over each cell
	for ( int i = 0 ; i < num_cells; i++)
	{
		cell_index = voronoi->index[i];
		num_cell_verticies = voronoi->num_cell_verticies[i];

		iv_zero(cell_sf_index);
		m_zero(bI);
		length_cell_index = 0;

		area = 0;

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
		double seg_length_m = sqrt(pow(v2[0]-v1[0],2) + pow(v2[1]-v1[1],2));  

		n_m[0] = (normalFactor)*1.00*(v2[1] - v1[1])/seg_length_m;
		n_m[1] = (normalFactor)*-1.00*(v2[0] - v1[0])/seg_length_m;

		// construct smoothed shape function gradients in the cell
		for ( int k = 0 ; k < num_cell_verticies ; k++)
		{
			// node 
			int indx;

			// create each segment v1 --- v2
			v1 = cell_verticies->me[cell_index[k]];
			if ( k < num_cell_verticies -1){
				v2 = cell_verticies->me[cell_index[k+1]];
				v3 = cell_verticies->me[cell_index[k+2]];		
				indx = k+1;
			}else{
				v2 = cell_verticies->me[cell_index[0]];
				v3 = cell_verticies->me[cell_index[1]];

				indx = 0;
			}

			// find length of segment 
			seg_length_m_1 = sqrt(pow(v3[0]-v2[0],2) + pow(v3[1]-v2[1],2));  

			// find outward normal of cell segment
			n_m_1[0] = (normalFactor)*1.00*(v3[1] - v2[1])/seg_length_m_1;
			n_m_1[1] = (normalFactor)*-1.00*(v3[0] - v2[0])/seg_length_m_1;


			area += (v1[0]*v2[1] - v2[0]*v1[1]);
			center[0] += v2[0];
			center[1] += v2[1];

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
					bI->me[idx][0] += 0.5*(n_m[0]*seg_length_m + n_m_1[0]*seg_length_m_1)*phi->ve[j];
					bI->me[idx][1] += 0.5*(n_m[1]*seg_length_m + n_m_1[1]*seg_length_m_1)*phi->ve[j];

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
					bI->me[length_cell_index-1][0] += 0.5*(n_m[0]*seg_length_m + n_m_1[0]*seg_length_m_1)*phi->ve[j];
					bI->me[length_cell_index-1][1] += 0.5*(n_m[1]*seg_length_m + n_m_1[1]*seg_length_m_1)*phi->ve[j];

				}
			}
			// int factor = 0.5 * ( n_im * l_m + n_i(m+1) * l_(m+1) ) 

			n_m[0] = n_m_1[0]; n_m[1] = n_m_1[1];
			seg_length_m = seg_length_m_1;
			center[0] = center[0]/num_cell_verticies;
			center[1] = center[1]/num_cell_verticies;

		} // end of loop over cell verticies


		area = 0.5*fabs(area);

		total_area += area;


	} // end of loop over cells

	//m_foutput(stdout, bI);
	//iv_foutput(stdout,cell_sf_index);

	printf("total area = %lf \n ",total_area);
	printf("Finished looping over cells and construction B \n ");

	// get shape function for all the verticies

	// if axi get shape function for all the nodes as well

	// get shape function for all the nodes


	


	/* Get list of points = verticies and node */

	// EFG_SF * _phi = malloc(1*sizeof(EFG_SF));


	// if (is_AXI == 1)
	// {

	// }
	// /*  Initialise area and length*/
	// scni->area = 0;
	// double seglength = 0;
	// scni->center = malloc(2*sizeof(double));
	// scni->center[0] = 0;
	// scni->center[1] = 0;
	// double normal[2] = {0,0};

	// /*  Recursive rule of integration given by Chen (2001) as
	//  *	bIi(x_L) = sum^Ns [ _phi_I(x_L^m)n_IL^m l_LM/2 + _phi_I(x_L^M+1) nIl^m lLM/2]
	//  *  */
	// /*  Start loop */
	// MAT * bI = m_get(efgBlock->numnode,2);

	// // check orientation of the polygon
	// gpc_vertex a = vorCell->contour->vertex[0];
	// gpc_vertex b = vorCell->contour->vertex[1];
	// gpc_vertex c = vorCell->contour->vertex[2];


	// int isCCW_ = isCCW(a,b,c);
	// int normalFactor;
	// // needed to make sure the normal is outward pointing when doing contour integration
	// if ( isCCW_ == 1)
	// {
	// 	normalFactor = 1;
	// }else{
	// 	normalFactor = -1;
	// }

	// gpc_vertex v1;
	// gpc_vertex v2;



	return 0;
}


int free_scni(SCNI * _scni){


	return 0;
}
