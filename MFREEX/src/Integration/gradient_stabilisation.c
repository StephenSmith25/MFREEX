#include "Integration/gradient_stabilisation.h"


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


stabalised_gradient **  gradient_stabilisation( shape_function_container * sf_point, meshfreeDomain * mfree )
{

	VEC * di = mfree->di;
	int num_points = mfree->num_nodes;
	int dim = mfree->dim;

	MAT * phi_der;
	IVEC * neighbours;
	MAT * phi_der_at_L;
	IVEC * neighbours_L;
	MAT * phi_der_der ;



	VEC * phi_der_cor = v_get(dim);
	VEC * g_l = v_get(dim);
	VEC * g_k = v_get(dim);
	VEC * g_kk = v_get(dim*dim);
	MAT * Hessian = m_get(dim,dim);
	MAT * a_Hessian_i = m_get(dim,dim);
	MAT * a_Hessian = m_get(dim,dim);

	int num_elements = dim*dim;
	int num_neighbours =0;

	stabalised_gradient ** G_s = malloc(num_points*sizeof(stabalised_gradient*));



	MAT * nodes = mfree->nodes;


	// loop over points
	for ( int i = 0 ; i < num_points ; i++)
	{ 
		G_s[i] = malloc(sizeof(stabalised_gradient));
		phi_der = sf_point->sf_list[i]->dphi;


		G_s[i]->g = m_get(phi_der->m,phi_der->n);

		phi_der_der = sf_point->sf_list[i]->d2phi;


		neighbours = sf_point->sf_list[i]->neighbours;
		num_neighbours = neighbours->max_dim;

		for ( int k = 0 ; k < num_neighbours ; k++){
			

			m_zero(a_Hessian);
			int indx_k = neighbours->ive[k];
			g_l = get_row(phi_der, k, g_l);
			g_kk = get_row(phi_der_der,k,g_kk);
			// smoothing length
			double h = di->ve[k];


			if ( dim == 1)
			{
				Hessian->me[0][0] = 0;
			}
			else if ( dim == 2)
			{
				Hessian->me[0][0] = g_kk->ve[0];					
				Hessian->me[0][1] = g_kk->ve[1];
				Hessian->me[1][0] = g_kk->ve[2];
				Hessian->me[1][1] = g_kk->ve[3];

			}else{
				Hessian->me[0][0] = g_kk->ve[0];					
				Hessian->me[0][1] = g_kk->ve[1];
				Hessian->me[0][2] = g_kk->ve[2];
				Hessian->me[1][0] = g_kk->ve[3];
				Hessian->me[1][1] = g_kk->ve[4];					
				Hessian->me[1][2] = g_kk->ve[5];
				Hessian->me[2][0] = g_kk->ve[6];
				Hessian->me[2][1] = g_kk->ve[7];
				Hessian->me[2][2] = g_kk->ve[7];

			}

			for ( int l = 0 ; l < num_neighbours ; l++)
			{
				// get neighbours at x_L
				int index_l = neighbours->ive[l];
				neighbours_L = sf_point->sf_list[index_l]->neighbours;
				// get phi derivatves at x_l
				phi_der_at_L = sf_point->sf_list[index_l]->dphi;
				int indx = findInt(indx_k,neighbours_L->ive,neighbours_L->max_dim);
				g_l = get_row(phi_der, l, g_l);
				//printf("indx = %d for k = %d and l = %d \n", indx, indx_k,index_l);
				if ( indx > -1){
					g_k = get_row(phi_der_at_L, indx, g_k);
					m_zero(a_Hessian_i);	
					v_outer_product(g_k, g_l, a_Hessian_i);
					//m_foutput(stdout, a_Hessian_i);
					__add__(a_Hessian->base, a_Hessian_i->base, a_Hessian->base, num_elements);
				}
			}

			// (H(phi_k) - grad(grad(phi)))
			__sub__(Hessian->base, a_Hessian->base, Hessian->base, num_elements);
			// (H(phi_k) - grad(grad(phi)))*h = g_s
			if ( dim == 1)
			{
				phi_der_cor->ve[0] = Hessian->me[0][0]*h;
			}else if ( dim == 2){
				phi_der_cor->ve[0] = (Hessian->me[0][0] + Hessian->me[0][1])*h;
				phi_der_cor->ve[1] = (Hessian->me[1][0] + Hessian->me[1][1])*h;
			}else{
				phi_der_cor->ve[0] = (Hessian->me[0][0] + Hessian->me[0][1] + Hessian->me[0][2])*h;
				phi_der_cor->ve[1] = (Hessian->me[1][0] + Hessian->me[1][1] + Hessian->me[1][2])*h;			
				phi_der_cor->ve[2] = (Hessian->me[2][0] + Hessian->me[2][1] + Hessian->me[2][2])*h;
			}

			set_row(G_s[i]->g, k, phi_der_cor);
		}
	}


	M_FREE(a_Hessian_i);
	M_FREE(a_Hessian);
	M_FREE(Hessian);
	V_FREE(g_kk);
	V_FREE(g_k);
	V_FREE(g_l);
	V_FREE(phi_der_cor);

	return G_s;
}
