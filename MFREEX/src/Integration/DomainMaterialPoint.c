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
#include "cellSearch.h"
#include "dsyev2.h"
#include "math/cmath3d.h"
#ifndef DIM 
#define DIM 2 
#endif 


// n choose 3, n is atleast 3
static int nChoose3[13] = {
	1, // base case n = 3 
	4, // n = 4;
	10,
	20,
	35,
	56,
	84,
	120,
	165,
	220,
	286,
	364,
	455,
};



void find_combinations(int combinations[][3],IVEC * points, int num_points)
{
    int i, j, k;
    int count = 0;
    for (i = 0; i < num_points - 2; i++)
    {
        for (j = i + 1; j < num_points - 1; j++)
        {
            for (k = j + 1; k < num_points; k++)
            {
                combinations[count][0] =  points->ive[i];
                combinations[count][1] =  points->ive[j];
                combinations[count][2] =  points->ive[k];
			
				++count;
			}
        }

    }

    return ;
}

static inline double  crossProduct(double * p1 , double * p2)
{
	return p1[0]*p2[1] - p1[1]*p2[0];
}

static inline bool SameSide(double p1[2], double p2[2], double a[2], double b[2] )
{
	double b_a[2] = {b[0] - a[0], b[1] - a[1] };
	double p1_a[2] = {p1[0] - a[0], p1[1] - a[1]};
	double p2_a[2] = {p2[0] - a[0], p2[1] - a[1]};

	double res_1 = crossProduct(b_a, p1_a);
	double res_2 = crossProduct(b_a, p2_a);
	double dot_res = res_1*res_2;


	if ( dot_res >= 0 )
		return true;
	else
		return false;
}

static inline bool PointInTriangle (double * p, double * a, double * b, double * c)
{


	if (SameSide(p,a, b,c) && SameSide(p,b, a,c) && SameSide(p,c, a,b) )
        return true;
    else
		return false;
}


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

int setDomainMaterialPoint(MAT * nodes, CELLS * cells, MATERIAL_POINT * MP)
{

	int num_nodes = nodes->m;
	int dim = nodes->n;



			VEC * distances = v_get(num_nodes); 
			PERM * order = px_get(num_nodes);
	


			// e.g find the 5 closest neighbours to x_p;
			for (int j = 0; j < num_nodes; ++j)
			{
				// find distance to point;
				distances->ve[j] = sq_distance(MP->coords_n_1, nodes->me[j], dim);


			}


			v_sort(distances,order);

	
			printf("distances = %lf \n", distances->ve[4]);

			MP->r_cutoff = MP->beta *  distances->ve[4];



			// Initially spherical domain of influence: 
			// MI = [r^2 0 0 ; 0 r^2 0 ; 0 0 r^2]
#if DIM == 2 
		MP->MI->me[0][0] =  MP->beta * MP->beta * distances->ve[4] *  distances->ve[4];
		MP->MI->me[1][1] =  MP->beta * MP->beta * distances->ve[4] *  distances->ve[4];
		MP->invMI->me[0][0] = 1.00/ MP->MI->me[0][0] ;
		MP->invMI->me[1][1] = 1.00/ MP->MI->me[1][1] ;

#elif DIM == 3
			MP->MI->me[0][0] = distances->ve[8];
			MP->MI->me[1][1] = distances->ve[8];
			MP->MI->me[2][2] = distances->ve[8];
			MP->invMI->me[0][0] = 1.00/distances->ve[8];
			MP->invMI->me[1][1] = 1.00/distances->ve[8];
			MP->invMI->me[2][2] = 1.00/distances->ve[8];
#endif 


			// Find neighbours 
			RangeSearchMaterialPoint(MP, nodes, cells);


			PX_FREE(order);
			V_FREE(distances);


			MP->l_0 = sqrt(MP->MI->me[0][0]);
			MP->l_1 = sqrt(MP->MI->me[1][1]);

	

			printf("GOT TO THE END OF THIS CALL\n");



	return 0;

}

int updateDomainMaterialPoint(MAT * nodes, CELLS * cells,  MATERIAL_POINT * MP)
{


	MAT * F = MP->stateNew->F;
	MAT * C = MP->stateNew->C;




	//mtrm_mlt(F,F,C);
	mmtr_mlt(F, F,C);


	// C->me[0][1] = 0;
	// C->me[1][0] = 0;

	// Find eigen values of network strain
	//m_inverse_small(MP->inc_F,MP->stateNew->m_temp4);
	double lambda1 = 0;
	double lambda2 = 0;
	double cs = 0;
	double sn = 0;

	dsyev2(C->me[0][0], C->me[0][1], C->me[1][1], &lambda1, &lambda2,
                   &cs,&sn);

	lambda1 = sqrt(lambda1);
	lambda2 = sqrt(lambda2);



	// if ( MP->index == 1412)
	// {
	// 	m_foutput(stdout,F );
	// 	m_foutput(stdout,C );

	// 	printf("cs and sin = %lf %lf \n",cs,sn);

	// 	exit(0);
	// }


	// using diagonal elements of C and ignoring
	// double l_0 = MP->l_0*C->me[0][0];
	// double l_1 = MP->l_1*C->me[1][1];

	// Find new lengths of the ellipse
	double l_0 = MP->l_0*lambda1;
	double l_1 = MP->l_1*lambda2;

	double MAX_ASPECT_RATIO = 3;
	double max_l0_l1 = max(l_0,l_1);
	double aspect_ratio = 0;
	if ( max_l0_l1 == l_0)
	{
		aspect_ratio = l_0/l_1;

		if ( aspect_ratio > MAX_ASPECT_RATIO)
		{
			l_1 = l_0/MAX_ASPECT_RATIO;
		}	

	}else{
		aspect_ratio = l_1/l_0;

		if ( aspect_ratio > MAX_ASPECT_RATIO)
			l_0 = l_1/MAX_ASPECT_RATIO;
	}

	

	// Find ellipse in global axis
	struct mat33 D = mdiag(1.00/(l_0*l_0),1.00/(l_1*l_1),0);
	struct mat33 R = mzero();

	R.m[0][0] = cs;
	R.m[0][1] = sn;
	R.m[1][0] = -sn;
	R.m[1][1] = cs;
	R.m[2][2] = 1;

	double theta = acos(cs);
	MP->theta = theta;

	struct mat33 RD = mtrmul(R, D);
	struct mat33 A = mmul(RD, R);

	MP->invMI->me[0][0] = A.m[0][0];
	MP->invMI->me[0][1] = A.m[0][1];
	MP->invMI->me[1][0] = A.m[1][0];
	MP->invMI->me[1][1] = A.m[1][1];



	RangeSearchMaterialPoint(MP, nodes, cells);



	return 0;

}


// // mtrm_mlt(F,MP->MI,MP->invMI);
	// // m_mlt(MP->invMI,F,MP->MI);
	// m_inverse_small(MP->MI,MP->invMI);



// 	double distance[MP->num_neighbours];

// 	// find min_num_neighbours for material point x_p
// 	for (j = 0; j < MP->num_neighbours; ++j)
// 	{
// 			// find distance to point;
// 			int index = MP->neighbours->ive[j];
// 			distance[j] = sq_distance(MP->coords_n_1, nodes->me[index], dim);

// 	}
// 	double a;
// 	for (i = 0 ; i < MP->num_neighbours ; i++)
// 	{
// 		for (  j = i+1 ; j < MP->num_neighbours ; j++)
// 		{
// 			if (distance[i] > distance[j])
// 			{
// 				a = distance[i];
// 				distance[i] = distance[j];
// 				distance[j] = a;
// 			}
// 		}
// 	}

// 	// update supports
// 	MP->r_cutoff = MP->beta *distance[6];

// #if DIM == 2 
// 	MP->MI->me[0][0] = MP->r_cutoff *  MP->r_cutoff ;
// 	MP->MI->me[1][1] =  MP->r_cutoff  *  MP->r_cutoff ;
// 	MP->invMI->me[0][0] = 1.00/ MP->MI->me[1][1] ;
// 	MP->invMI->me[1][1] = 1.00/ MP->MI->me[1][1] ;

// #elif DIM == 3
// 	MP->MI->me[0][0] = distance[8];
// 	MP->MI->me[1][1] = distance[8];
// 	MP->MI->me[2][2] = distance[8];
// 	MP->invMI->me[0][0] = 1.00/distance[8];
// 	MP->invMI->me[1][1] = 1.00/distance[8];
// 	MP->invMI->me[2][2] = 1.00/distance[8];
// #endif 



	// update neighbours 
	// MP->num_neighbours = neighbour_RangeSearch(MP->neighbours
	//  		,cells, MP->coords_n_1, MP->r_cutoff, nodes);