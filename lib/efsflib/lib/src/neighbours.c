#include "neighbours.h"

#include "matrix.h"
#include "matrix2.h"
#include <stdbool.h>


const double TOLEREANCE = 1e-12;

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


static inline double rectangle_area(double x1, double y1, double x2, double y2, 
                            double x3, double y3) 
{ 
    return fabs((x1 * (y2 - y3) + x2 * (y3 - y1) +  
                x3 * (y1 - y2)) / 2.0); 
} 


static inline bool inside_rectangle(double x1, double y1, double x2, double y2, double x3, double y3,
	double x4, double y4, double x, double y)
{
	
  /* Calculate area of rectangle ABCD */
    double A = rectangle_area(x1, y1, x2, y2, x3, y3) +  
              rectangle_area(x1, y1, x4, y4, x3, y3); 
  
    /* Calculate area of triangle PAB */
    double A1 = rectangle_area(x, y, x1, y1, x2, y2); 
  
    /* Calculate area of triangle PBC */
    double A2 = rectangle_area(x, y, x2, y2, x3, y3); 
  
    /* Calculate area of triangle PCD */
    double A3 = rectangle_area(x, y, x3, y3, x4, y4); 
  
    /* Calculate area of triangle PAD */
    double A4 = rectangle_area(x, y, x1, y1, x4, y4); 
  
    /* Check if sum of A1, A2, A3 and A4  
       is same as A */
    double A_sum = A1+A2+A3+A4;

    if ( ((fabs(A-A_sum))) < TOLEREANCE)
    {
    	return true;

    }else{
    	return false;
    }

}


IVEC * get_point_neighbours(IVEC * neighbours, double * x, meshfreeDomain * mfree)
{

	int dim = mfree->dim;
	int numnodes = mfree->num_nodes;
	MAT * nodes = mfree->nodes;



	switch(mfree->kernel_support)
	{
		case (RADIAL):
		{


		VEC * domainSize = mfree->di;

	
	
		double distance = 0;
		int count_neighbours = 0;


		for (int i = 0; i < numnodes; ++i)
		{
			distance = sq_distance(x, nodes->me[i], dim);


			if ( distance < domainSize->ve[i])
			{
				neighbours->ive[count_neighbours] = i;
				count_neighbours = count_neighbours + 1;

				if ( count_neighbours >= neighbours->max_dim)
				{

					IVEC * temp = iv_get(neighbours->max_dim);
					iv_copy(neighbours,temp);
					iv_resize(neighbours,neighbours->max_dim+5);

					for (int k = 0 ; k < temp->max_dim; ++k)
					{
						neighbours->ive[k] = temp->ive[k];
					}
					IV_FREE(temp);

				}
			}
		}
		if ( neighbours->max_dim > count_neighbours)
		{
			IVEC * temp = iv_get(neighbours->max_dim);
			iv_copy(neighbours,temp);
			iv_resize(neighbours,count_neighbours);
			neighbours->max_dim = count_neighbours;

			for (int k = 0 ; k < count_neighbours; ++k)
			{
				neighbours->ive[k] = temp->ive[k];
			}
			IV_FREE(temp);

		}


		break;


	}
	case (RECTANGULAR):
	{

		bool is_inside; 
		int count_neighbours = 0;

		double x1,x2,x3,x4;
		double y1,y2,y3,y4;
	
		double dmx,dmy,dmz;
		double x_q,y_q,z_q;
		double x_I,y_I,z_I;
		x_q = x[0];
		y_q = x[1];
		for (int i = 0; i < numnodes; ++i)
		{

			// Form x1 x2 x3 x4 of each rectangle
			

			if ( dim == 2)
			{
				x_I = nodes->me[i][0];
				y_I = nodes->me[i][1];
				dmx = mfree->di_tensor->me[i][0];
				dmy = mfree->di_tensor->me[i][1];



				x1 = x_I - dmx;
				y1 = y_I-dmy;

				x2 = x_I + dmx;
				y2 = y_I - dmy;

				x3 = x_I+dmx;
				y3 = y_I+dmy;

				x4 = x_I - dmx;
				y4 = y_I+dmy;


				is_inside = inside_rectangle(x1, y1,x2,y2, x3, y3,x4, y4, x_q,y_q);


			}else if ( dim ==3)
			{

				fprintf(stderr,"need to implement this\n");

			}


			if ( is_inside == true)
			{
				neighbours->ive[count_neighbours] = i;
				count_neighbours = count_neighbours + 1;

				if ( count_neighbours >= neighbours->max_dim)
				{

					IVEC * temp = iv_get(neighbours->max_dim);
					iv_copy(neighbours,temp);
					iv_resize(neighbours,neighbours->max_dim+5);

					for (int k = 0 ; k < temp->max_dim; ++k)
					{
						neighbours->ive[k] = temp->ive[k];
					}
					IV_FREE(temp);

				}
			}
		}
		if ( neighbours->max_dim > count_neighbours)
		{
			IVEC * temp = iv_get(neighbours->max_dim);
			iv_copy(neighbours,temp);
			iv_resize(neighbours,count_neighbours);
			neighbours->max_dim = count_neighbours;

			for (int k = 0 ; k < count_neighbours; ++k)
			{
				neighbours->ive[k] = temp->ive[k];
			}
			IV_FREE(temp);

		}


		break;
	}
	case (ELLIPTICAL):
	{

		double distance = 0;
		int count_neighbours = 0;
		double beta = mfree->beta;
		double d_I_s = 0;
		MAT * MI ;


		double M11 = 0;
		double M12 = 0;
		double M21 = 0;
		double M22 = 0;



		for (int i = 0; i < numnodes; ++i)
		{
			distance = sq_distance(x, nodes->me[i], dim);
			double xS[2] = {x[0] - nodes->me[i][0], x[1] - nodes->me[i][1]};
			MI = mfree->MI[i];


			M11 = MI->me[0][0];
			M12 = MI->me[0][1];
			M21 = MI->me[1][0];
			M22 = MI->me[1][1];
			distance = xS[0]*(M11*xS[0] + M12 * xS[1]) + xS[1]*(M21*xS[0] + M22 * xS[1]);



			if ( distance <= 1)
			{
				neighbours->ive[count_neighbours] = i;
				count_neighbours = count_neighbours + 1;

				if ( count_neighbours >= neighbours->max_dim)
				{

					IVEC * temp = iv_get(neighbours->max_dim);
					iv_copy(neighbours,temp);
					iv_resize(neighbours,neighbours->max_dim+5);

					for (int k = 0 ; k < temp->max_dim; ++k)
					{
						neighbours->ive[k] = temp->ive[k];
					}
					IV_FREE(temp);

				}
			}
		}
		if ( neighbours->max_dim > count_neighbours)
		{
			IVEC * temp = iv_get(neighbours->max_dim);
			iv_copy(neighbours,temp);
			iv_resize(neighbours,count_neighbours);
			neighbours->max_dim = count_neighbours;

			for (int k = 0 ; k < count_neighbours; ++k)
			{
				neighbours->ive[k] = temp->ive[k];
			}
			IV_FREE(temp);

		}

		break;
	}
	}

	return neighbours;

}



IVEC * point_neighbours(double * x, meshfreeDomain * mfree)
{

	int dim = mfree->dim;
	int numnodes = mfree->num_nodes;
	IVEC * neighbours = iv_get((int)(floor ( sqrt(numnodes))));
	MAT * nodes = mfree->nodes;



	switch(mfree->kernel_support)
	{
		case (RADIAL):
		{


		VEC * domainSize = mfree->di;

	
	
		double distance = 0;
		int count_neighbours = 0;


		for (int i = 0; i < numnodes; ++i)
		{
			distance = sq_distance(x, nodes->me[i], dim);


			if ( distance < domainSize->ve[i])
			{
				neighbours->ive[count_neighbours] = i;
				count_neighbours = count_neighbours + 1;

				if ( count_neighbours >= neighbours->max_dim)
				{

					IVEC * temp = iv_get(neighbours->max_dim);
					iv_copy(neighbours,temp);
					iv_resize(neighbours,neighbours->max_dim+5);

					for (int k = 0 ; k < temp->max_dim; ++k)
					{
						neighbours->ive[k] = temp->ive[k];
					}
					IV_FREE(temp);

				}
			}
		}
		if ( neighbours->max_dim > count_neighbours)
		{
			IVEC * temp = iv_get(neighbours->max_dim);
			iv_copy(neighbours,temp);
			iv_resize(neighbours,count_neighbours);
			neighbours->max_dim = count_neighbours;

			for (int k = 0 ; k < count_neighbours; ++k)
			{
				neighbours->ive[k] = temp->ive[k];
			}
			IV_FREE(temp);

		}


		break;


	}
	case (RECTANGULAR):
	{

		bool is_inside; 
		int count_neighbours = 0;

		double x1,x2,x3,x4;
		double y1,y2,y3,y4;
	
		double dmx,dmy,dmz;
		double x_q,y_q,z_q;
		double x_I,y_I,z_I;
		x_q = x[0];
		y_q = x[1];
		for (int i = 0; i < numnodes; ++i)
		{

			// Form x1 x2 x3 x4 of each rectangle
			

			if ( dim == 2)
			{
				x_I = nodes->me[i][0];
				y_I = nodes->me[i][1];
				dmx = mfree->di_tensor->me[i][0];
				dmy = mfree->di_tensor->me[i][1];



				x1 = x_I - dmx;
				y1 = y_I-dmy;

				x2 = x_I + dmx;
				y2 = y_I - dmy;

				x3 = x_I+dmx;
				y3 = y_I+dmy;

				x4 = x_I - dmx;
				y4 = y_I+dmy;


				is_inside = inside_rectangle(x1, y1,x2,y2, x3, y3,x4, y4, x_q,y_q);


			}else if ( dim ==3)
			{

				fprintf(stderr,"need to implement this\n");

			}


			if ( is_inside == true)
			{
				neighbours->ive[count_neighbours] = i;
				count_neighbours = count_neighbours + 1;

				if ( count_neighbours >= neighbours->max_dim)
				{

					IVEC * temp = iv_get(neighbours->max_dim);
					iv_copy(neighbours,temp);
					iv_resize(neighbours,neighbours->max_dim+5);

					for (int k = 0 ; k < temp->max_dim; ++k)
					{
						neighbours->ive[k] = temp->ive[k];
					}
					IV_FREE(temp);

				}
			}
		}
		if ( neighbours->max_dim > count_neighbours)
		{
			IVEC * temp = iv_get(neighbours->max_dim);
			iv_copy(neighbours,temp);
			iv_resize(neighbours,count_neighbours);
			neighbours->max_dim = count_neighbours;

			for (int k = 0 ; k < count_neighbours; ++k)
			{
				neighbours->ive[k] = temp->ive[k];
			}
			IV_FREE(temp);

		}


		break;
	}
	case (ELLIPTICAL):
	{

		double distance = 0;
		int count_neighbours = 0;
		double beta = mfree->beta;
		double d_I_s = 0;
		MAT * MI ;


		double M11 = 0;
		double M12 = 0;
		double M21 = 0;
		double M22 = 0;



		for (int i = 0; i < numnodes; ++i)
		{
			distance = sq_distance(x, nodes->me[i], dim);
			double xS[2] = {x[0] - nodes->me[i][0], x[1] - nodes->me[i][1]};
			MI = mfree->MI[i];


			M11 = MI->me[0][0];
			M12 = MI->me[0][1];
			M21 = MI->me[1][0];
			M22 = MI->me[1][1];
			distance = xS[0]*(M11*xS[0] + M12 * xS[1]) + xS[1]*(M21*xS[0] + M22 * xS[1]);



			if ( distance <= 1)
			{
				neighbours->ive[count_neighbours] = i;
				count_neighbours = count_neighbours + 1;

				if ( count_neighbours >= neighbours->max_dim)
				{

					IVEC * temp = iv_get(neighbours->max_dim);
					iv_copy(neighbours,temp);
					iv_resize(neighbours,neighbours->max_dim+5);

					for (int k = 0 ; k < temp->max_dim; ++k)
					{
						neighbours->ive[k] = temp->ive[k];
					}
					IV_FREE(temp);

				}
			}
		}
		if ( neighbours->max_dim > count_neighbours)
		{
			IVEC * temp = iv_get(neighbours->max_dim);
			iv_copy(neighbours,temp);
			iv_resize(neighbours,count_neighbours);
			neighbours->max_dim = count_neighbours;

			for (int k = 0 ; k < count_neighbours; ++k)
			{
				neighbours->ive[k] = temp->ive[k];
			}
			IV_FREE(temp);

		}

		break;
	}
	}

	return neighbours;

}
