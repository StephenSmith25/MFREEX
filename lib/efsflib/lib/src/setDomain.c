#include "setDomain.h"




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

			double dmax = mfree->dmax_radial;
			double distance = 0;
			double distance_min = 1e6;
			// all nodes have same support size

			for (int i = 0; i < num_nodes; ++i)
			{
			distance_min = 1e6;
			distance = 1e6;


			// find distance from node I to each node J 
			for (int j = 0; j < num_nodes; ++j)
			{
				if ( i != j){
					// find distance to point;
					distance = sq_distance(nodes->me[i], nodes->me[j], dim);
					if (distance < distance_min )
					{
						distance_min = distance;
					}
				}
			}
			mfree->di->ve[i] = distance_min*dmax;
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

		int min_num_neighbours = 6;


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


			break;
		}
		case(ELLIPTICAL):
		{
			break;
		}
		default:

		fprintf(stderr,"kernel shape not set\n");
	}


	return 0;

}


