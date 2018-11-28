#include "v_outer_product.h"




MAT * v_outer_product(VEC * v1, VEC * v2, MAT * out)
{

	if ( v1->max_dim != v2->max_dim)
	{
		fprintf(stderr,"outer product not valid for vectors of different sizes \n");
	}

	int dim = v1->max_dim;



	if (out == MNULL)
	{
		out = m_get(dim,dim);
	}else{
		if ( ( out->m != dim) || ( out->n != dim)){
			m_resize(out, v1->max_dim, v2->max_dim);
		}
	}

	for ( int i = 0 ; i < dim ; i++)
	{
		for (int j = 0; j < dim; ++j)
		{
			out->me[i][j] = v1->ve[i] * v2->ve[j];
		}

	}

	return out;
}
