#include "Boundary/Traction/Pressure/pressure_load.h"




pressure_boundary * new_pressure_boundary(IVEC * points,  meshfreeDomain * mFree)
{

	pressure_boundary * pB = malloc(1*sizeof(pressure_boundary));
	MAT * coords = mFree->nodes;
	int dim = coords->n;
	int num_points = points->max_dim;
	

	// points
	pB->points = points;
	pB->coords = m_get(num_points,dim);

	// set up traction coordinates
	for ( int i = 0 ; i < num_points ; i++)
	{
		pB->coords->me[i] = coords->me[pB->points->ive[i]];
	}





	// integration will be performed using a trapzoidal rule
	pB->sf_traction = mls_shapefunction(pB->coords, 1, mFree);

	pB->segment_normals = m_get(num_points-1,dim);
	pB->segment_weights = malloc((num_points-1)*sizeof(double));


	// find segment weights and normals
	double x1 = pB->coords->me[0][0], y1 = pB->coords->me[0][1];
	double x2, y2;
	double seg_length;
	for ( int i = 0 ; i < num_points -1  ; i++)
	{

		x2 = pB->coords->me[i+1][0], y2 = pB->coords->me[i+1][1];


		// length of segment
		seg_length =  sqrt(pow(x2-x1,2) + pow(y2-y1,2));


		if ( pB->is_axi == 1)
		{
			pB->segment_weights[i] = seg_length*2*M_PI*((x1+x2)/2.00);
		}else{
			pB->segment_weights[i] = seg_length;
		}

		// normal
		pB->segment_normals->me[i][0] = (y2-y1)/seg_length;
		pB->segment_normals->me[i][1] = -(x2-x1)/seg_length;

		x1 = x2;
		y1 = y2;

	}



	return pB;
}


int update_pressure_boundary(pressure_boundary *pB, MAT * coords)
{

	int num_points = pB->points->max_dim;

	// set up traction coordinates
	for ( int i = 0 ; i < num_points ; i++)
	{
		pB->coords->me[i] = coords->me[pB->points->ive[i]];
	}


	// find segment weights and normals
	double x1 = pB->coords->me[0][0], y1 = pB->coords->me[0][1];
	double x2, y2;
	double seg_length;
	for ( int i = 0 ; i < num_points -1  ; i++)
	{

		x2 = pB->coords->me[i+1][0], y2 = pB->coords->me[i+1][1];


		// length of segment
		seg_length =  sqrt(pow(x2-x1,2) + pow(y2-y1,2));



		if ( pB->is_axi == 1)
		{
			pB->segment_weights[i] = seg_length*2*M_PI*((x1+x2)/2.00);
		}else{
			pB->segment_weights[i] = seg_length;
			
		}


		// normal
		pB->segment_normals->me[i][0] = (y2-y1)/seg_length;
		pB->segment_normals->me[i][1] = -(x2-x1)/seg_length;

		x1 = x2;
		y1 = y2;

	}
	return 0;
}

int assemble_pressure_load(VEC * Fext, double Pressure, pressure_boundary * pB)
{
	double weight = 0;	
	int numPoints = pB->points->max_dim;
	double surface_traction[2];
	double * segment_normal;

	IVEC * neighbours_n1;
	VEC * phi_n1;

	IVEC * neighbours_n2;
	VEC * phi_n2;

	int num_neighbours_n1;
	int num_neighbours_n2;

	double intFactor = 0;
	//double area_sum = 0;

	for ( int i = 0 ; i < numPoints - 1; i++ )
	{			

			phi_n1 = pB->sf_traction->sf_list[i]->phi;
			neighbours_n1 = pB->sf_traction->sf_list[i]->neighbours;
			num_neighbours_n1 = neighbours_n1->max_dim;
			phi_n2 = pB->sf_traction->sf_list[i+1]->phi;
			neighbours_n2 = pB->sf_traction->sf_list[i+1]->neighbours;
			num_neighbours_n2 = neighbours_n2->max_dim;

			segment_normal = pB->segment_normals->me[i];
			// x component of traction
			surface_traction[0] = Pressure*segment_normal[0];
			// y component of traction
			surface_traction[1] = Pressure*segment_normal[1];
			int MAX_NEIGHBOURS = max(num_neighbours_n1,num_neighbours_n2);
			intFactor = pB->segment_weights[i];

			for ( int k = 0 ; k < MAX_NEIGHBOURS; k++)
			{

				if ( k < num_neighbours_n1)
				{
					Fext->ve[2*neighbours_n1->ive[k]] += 0.5*intFactor*surface_traction[0]*phi_n1->ve[k] ;
					Fext->ve[2*neighbours_n1->ive[k]+1] += 0.5*intFactor*surface_traction[1]*phi_n1->ve[k]  ;

				}
				if ( k < num_neighbours_n2)
				{
					Fext->ve[2*neighbours_n2->ive[k]] += 0.5*intFactor*surface_traction[0]*phi_n2->ve[k] ;
					Fext->ve[2*neighbours_n2->ive[k]+1] += 0.5*intFactor*surface_traction[1]*phi_n2->ve[k] ;

				}

			}


			//area_sum +=  pB->segment_weights[i];

	




	}

	// internal area 
	//printf("internal area = %lf \n", area_sum);


	return 0;
}

int free_pressure_boundary(pressure_boundary * pB)
{
	return 0;
}
