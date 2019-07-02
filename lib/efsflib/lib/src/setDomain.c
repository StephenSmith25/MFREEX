#include "setDomain.h"


const int min_num_neighbours = 4;

static inline double sq_distance(double *point_1, double * point_2, int  dim)
{
	double distance = 0;
	for (int i = 0; i < dim; ++i)
	{
		distance += pow(point_2[i] - point_1[i], 2);
		/* code */
	}

	return sqrt(distance);

}



int setDomain(meshfreeDomain * mfree)
{

	MAT * nodes = mfree->nodes;
	int constant_support_size = mfree->is_constant_support_size;

	int num_nodes = mfree->num_nodes;
	int dim = mfree->dim;




	if ( strcmp(mfree->kernel_shape,"rectangular") == 0)
	{
		mfree->kernel_support = RECTANGULAR;

	}else if (  strcmp(mfree->kernel_shape,"radial") == 0)
	{
		mfree->kernel_support = RADIAL;

	}else if (  strcmp(mfree->kernel_shape,"elliptical") == 0)

		mfree->kernel_support = ELLIPTICAL;

	else{
	}



	switch(mfree->kernel_support)
	{
		case RADIAL:
		{


			VEC * distances = v_get(num_nodes); 
			PERM * order = px_get(num_nodes);
			double dmax = mfree->dmax_radial;
			// all nodes have same support size

			if ( mfree->di == VNULL)
			{
			mfree->di = v_get(num_nodes);
			}



			for (int i = 0; i < num_nodes; ++i)
			{
	



			// find distance from node I to each node J 
			for (int j = 0; j < num_nodes; ++j)
			{
				// find distance to point;
				distances->ve[j] = sq_distance(nodes->me[i], nodes->me[j], dim);

			}



			v_sort(distances,order);
			mfree->di->ve[i] = distances->ve[min_num_neighbours]*dmax;



			}	

			int index = 0;
			if  ( constant_support_size == 1)
			{

			double max_distance = v_max(mfree->di,&index);

			for (int i = 0; i < num_nodes; ++i)
			{
				mfree->di->ve[i] = max_distance;
				}
			}

			break;
		}
		case(RECTANGULAR):
		{


		if ( mfree->di_tensor == MNULL)
		{
		mfree->di_tensor = m_get(num_nodes,dim);
		}



		VEC * distances = v_get(num_nodes); 
		PERM * order = px_get(num_nodes);

		// Find min_num_neighbours of point I
		double di_k[dim];



		for (int i = 0; i < num_nodes; ++i)
			{


			// find distance from node I to each node J 
			for (int j = 0; j < num_nodes; ++j)
			{
				// find distance to point;
				distances->ve[j] = sq_distance(nodes->me[i], nodes->me[j], dim);

			}



			v_sort(distances,order);
			double di_k[3];
			double max_distance[dim];

			// find dmx dmz based on this

			for ( int j = 0 ; j < min_num_neighbours ; j++)
			{
				for ( int k = 0 ; k < dim ; k++){
				double max_distance_j = fabs(nodes->me[i][k] - nodes->me[order->pe[j]][k]);
				if ( max_distance_j > max_distance[k])
				{
					max_distance[k] = max_distance_j;
				}
				}
			}
			for ( int k = 0 ; k < dim ; k++)
			{

			 		mfree->di_tensor->me[i][k] = mfree->dmax_tensor[k]*max_distance[k];
			}




			}// end loop over nodes	




			if ( constant_support_size == 1)
			{
				double di_max[dim];
				int index = 0;

				for ( int k = 0 ; k < dim ; k++)
				{
					distances = get_col(mfree->di_tensor,k,distances);
					di_max[k] = v_max(distances,&index);

				}

				for ( int i = 0  ; i < num_nodes ; i++)
				{
					for ( int k = 0 ; k < dim ; k++)
					{
						mfree->di_tensor->me[i][k] = di_max[k];
					}
				}



		}

			PX_FREE(order);
			V_FREE(distances);





			break;
		}
		case(ELLIPTICAL):
		{

			if ( mfree->MI == NULL)
			{
				mfree->MI = malloc(mfree->num_nodes*sizeof(MAT *));
				for (int i = 0; i < num_nodes; ++i)
				{
				
				mfree->MI[i] = m_get(2,2);
				}
			}


			VEC * distances = v_get(num_nodes); 
			PERM * order = px_get(num_nodes);

			// form moment matrix

			if  ( dim !=2)
			{
				fprintf(stderr,"ERROR: Eliptical basis only usable in 2D currently");
				return -1;
			}


			MAT * MI = m_get(2,2);

			for (int i = 0; i < num_nodes; ++i)
			{
				
				m_zero(MI);


				// find distance from node I to each node J 
				for (int j = 0; j < num_nodes; ++j)
				{
					// find distance to point;
					distances->ve[j] = sq_distance(nodes->me[i], nodes->me[j], dim);


				}

				// Form MI with the N neighbours
				v_sort(distances,order);
				double xI = mfree->nodes->me[i][0];
				double yI = mfree->nodes->me[i][1];
				double xJ,yJ;

				int k = 0;

				double distance_n = -1;
				double distance_n_1 = -1;

				int j = 0;
					while (( j < 100) && ( k < min_num_neighbours))
					{
						distance_n_1 = distances->ve[j];

						xJ = mfree->nodes->me[order->pe[j]][0];
						yJ = mfree->nodes->me[order->pe[j]][1];




						MI->me[0][0] += (xJ-xI)*(xJ-xI);
						MI->me[0][1] += (xJ-xI)*(yJ-yI);
						MI->me[1][0] += (xJ-xI)*(yJ-yI);
						MI->me[1][1] += (yJ-yI)*(yJ-yI);


						j++;
						




					if ( fabs(distance_n - distance_n_1) < 1e-9 )
					{

					}else{
						k = k+1;
					}

					// update counters
					distance_n = distance_n_1;
					j++;

					}


	
				sm_mlt(mfree->beta,MI,MI);
				mfree->MI[i] = m_inverse(MI,mfree->MI[i]);









			}


			m_free(MI);

			PX_FREE(order);
			V_FREE(distances);

			break;


		}
		default:

		printf("domain type not set\n");



	}


	return 0;

}



