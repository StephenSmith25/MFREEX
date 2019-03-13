/* ********************************
 * Author:       Stephen Smith 
 * License:	     Not yet decided 
 * Description:  Sets domain of influence for each material point
 *
 *//** @file DomainMaterialPoint.h *//*
 *
 ********************************/

#include "Integration/DomainMaterialPoint.h"
#include "ShapeFunction/neighbours_materialpoint.h"

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

int setDomainMaterialPoint(MAT * nodes, MATERIAL_POINT * MP)
{

	int num_nodes = nodes->m;
	int dim = nodes->n;
	int min_num_neighbours = 6;



	switch(MP->kernel_support)
	{
		case(RADIAL):
		{
	
			MP->MI = m_get(dim,dim);
			MP->invMI = m_get(dim,dim);


			VEC * distances = v_get(num_nodes); 
			PERM * order = px_get(num_nodes);
	

			// find min_num_neighbours for material point x_p

			// e.g find the 5 closest neighbours to x_p;
			for (int j = 0; j < num_nodes; ++j)
			{
				// find distance to point;
				distances->ve[j] = sq_distance(MP->coords_n_1, nodes->me[j], dim);


			}

			// Form MI with the N neighbours
			v_sort(distances,order);

			for ( int k = 0 ; k < dim ; k++)
			{
				MP->MI->me[k][k] = distances->ve[min_num_neighbours-1];
			}
			
	
			sm_mlt(MP->beta,MP->MI,MP->MI);



			MP->invMI = m_inverse_small(MP->MI,MP->invMI);


			PX_FREE(order);
			V_FREE(distances);

			break;
		}

		
		case(ELLIPTICAL):
		{

			// if ( mfree->MI == NULL)
			// {
			// 	mfree->MI = malloc(mfree->num_nodes*sizeof(MAT *));
			// 	for (int i = 0; i < num_nodes; ++i)
			// 	{
				
			// 	mfree->MI[i] = m_get(2,2);
			// 	}
			// }


			// VEC * distances = v_get(num_nodes); 
			// PERM * order = px_get(num_nodes);

			// // form moment matrix

			// if  ( dim !=2)
			// {
			// 	fprintf(stderr,"ERROR: Eliptical basis only usable in 2D currently");
			// 	return -1;
			// }


			// MAT * MI = m_get(2,2);

			// for (int i = 0; i < num_nodes; ++i)
			// {
				
			// 	m_zero(MI);


			// 	// find distance from node I to each node J 
			// 	for (int j = 0; j < num_nodes; ++j)
			// 	{
			// 		// find distance to point;
			// 		distances->ve[j] = sq_distance(nodes->me[i], nodes->me[j], dim);


			// 	}

			// 	// Form MI with the N neighbours
			// 	v_sort(distances,order);
			// 	double xI = mfree->nodes->me[i][0];
			// 	double yI = mfree->nodes->me[i][1];
			// 	double xJ,yJ;

			// 	int k = 0;

			// 	double distance_n = -1;
			// 	double distance_n_1 = -1;

			// 	int j = 0;
			// 		while (( j < 100) && ( k < min_num_neighbours))
			// 		{
			// 			distance_n_1 = distances->ve[j];

			// 			xJ = mfree->nodes->me[order->pe[j]][0];
			// 			yJ = mfree->nodes->me[order->pe[j]][1];




			// 			MI->me[0][0] += (xJ-xI)*(xJ-xI);
			// 			MI->me[0][1] += (xJ-xI)*(yJ-yI);
			// 			MI->me[1][0] += (xJ-xI)*(yJ-yI);
			// 			MI->me[1][1] += (yJ-yI)*(yJ-yI);


			// 			j++;
						




			// 		if ( fabs(distance_n - distance_n_1) < 1e-9 )
			// 		{

			// 		}else{
			// 			k = k+1;
			// 		}

			// 		// update counters
			// 		distance_n = distance_n_1;
			// 		j++;

			// 		}


	
			// 	sm_mlt(mfree->beta,MI,MI);
			// 	mfree->MI[i] = m_inverse(MI,mfree->MI[i]);









			// }


			// m_free(MI);

			// PX_FREE(order);
			// V_FREE(distances);

			// break;


		}
		default:

		printf("domain type not set\n");


	}


	return 0;

}

int updateDomainMaterialPoint(MAT * nodes, MATERIAL_POINT * MP)
{

	int num_nodes = nodes->m;
	int dim = nodes->n;
	int min_num_neighbours = 6;



	switch(MP->kernel_support)
	{
		case(RADIAL):
		{
	


			VEC * distances = v_get(num_nodes); 
			PERM * order = px_get(num_nodes);
	

			// find min_num_neighbours for material point x_p

			// e.g find the 5 closest neighbours to x_p;
			for (int j = 0; j < num_nodes; ++j)
			{
				// find distance to point;
				distances->ve[j] = sq_distance(MP->coords_n_1, nodes->me[j], dim);


			}

			// Form MI with the N neighbours
			v_sort(distances,order);

			m_zero(MP->MI);

			for ( int k = 0 ; k < dim ; k++)
			{
				MP->MI->me[k][k] = distances->ve[min_num_neighbours-1];
			}
			
	
			sm_mlt(MP->beta,MP->MI,MP->MI);



			MP->invMI = m_inverse_small(MP->MI,MP->invMI);
			
			MP->shape_function->neighbours = get_materialpoint_neighbours
			(MP->shape_function->neighbours,MP,nodes);

			PX_FREE(order);
			V_FREE(distances);

			break;
		}

		
		case(ELLIPTICAL):
		{

			// if ( mfree->MI == NULL)
			// {
			// 	mfree->MI = malloc(mfree->num_nodes*sizeof(MAT *));
			// 	for (int i = 0; i < num_nodes; ++i)
			// 	{
				
			// 	mfree->MI[i] = m_get(2,2);
			// 	}
			// }


			// VEC * distances = v_get(num_nodes); 
			// PERM * order = px_get(num_nodes);

			// // form moment matrix

			// if  ( dim !=2)
			// {
			// 	fprintf(stderr,"ERROR: Eliptical basis only usable in 2D currently");
			// 	return -1;
			// }


			// MAT * MI = m_get(2,2);

			// for (int i = 0; i < num_nodes; ++i)
			// {
				
			// 	m_zero(MI);


			// 	// find distance from node I to each node J 
			// 	for (int j = 0; j < num_nodes; ++j)
			// 	{
			// 		// find distance to point;
			// 		distances->ve[j] = sq_distance(nodes->me[i], nodes->me[j], dim);


			// 	}

			// 	// Form MI with the N neighbours
			// 	v_sort(distances,order);
			// 	double xI = mfree->nodes->me[i][0];
			// 	double yI = mfree->nodes->me[i][1];
			// 	double xJ,yJ;

			// 	int k = 0;

			// 	double distance_n = -1;
			// 	double distance_n_1 = -1;

			// 	int j = 0;
			// 		while (( j < 100) && ( k < min_num_neighbours))
			// 		{
			// 			distance_n_1 = distances->ve[j];

			// 			xJ = mfree->nodes->me[order->pe[j]][0];
			// 			yJ = mfree->nodes->me[order->pe[j]][1];




			// 			MI->me[0][0] += (xJ-xI)*(xJ-xI);
			// 			MI->me[0][1] += (xJ-xI)*(yJ-yI);
			// 			MI->me[1][0] += (xJ-xI)*(yJ-yI);
			// 			MI->me[1][1] += (yJ-yI)*(yJ-yI);


			// 			j++;
						




			// 		if ( fabs(distance_n - distance_n_1) < 1e-9 )
			// 		{

			// 		}else{
			// 			k = k+1;
			// 		}

			// 		// update counters
			// 		distance_n = distance_n_1;
			// 		j++;

			// 		}


	
			// 	sm_mlt(mfree->beta,MI,MI);
			// 	mfree->MI[i] = m_inverse(MI,mfree->MI[i]);









			// }


			// m_free(MI);

			// PX_FREE(order);
			// V_FREE(distances);

			// break;


		}
		default:

		printf("domain type not set\n");


	}


	return 0;

}