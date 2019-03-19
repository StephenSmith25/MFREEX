#include "ShapeFunction/mls_shapefunction_materialpoint.h"

shape_function * new_shape_function(int compute, char * basis_type, int dim)
{




	// get size of matricies based on basis used
	int dim_p;
	if ( strcmp(basis_type,"linear") == 0 )
	{
		if ( dim == 1){
			dim_p = 2;
		}else if ( dim == 2){
			dim_p = 3;
		}else{
			dim_p = 4;
		}
	}else if ( strcmp(basis_type,"quadratic") == 0)
	{
		if ( dim == 1){
			dim_p = 3;
		}else if ( dim == 2){
			dim_p = 6;
		}else{
			dim_p = 10;
		}

	}else {
		fprintf(stderr,"basis type not specified, or unknown");
	}	


	shape_function * sf_point = malloc(1*sizeof(shape_function));

		// 
	sf_point->phi = v_get(15);
	sf_point->dphi = m_get(10,10);


	// INTERNAL MATRICIES
	sf_point->A = m_get(dim_p,dim_p);
	sf_point->B = m_get(dim_p,dim_p);
	sf_point->basis_xi = m_get(dim_p,1);
	sf_point->ppT = m_get(dim_p,dim_p);
	sf_point->bi = v_get(dim_p);
	sf_point->p_xi = v_get(dim_p);
	sf_point->weights = v_get(1+dim);
	sf_point->basis = m_get(1,1);
	sf_point->p = v_get(dim_p);
	sf_point->invA = m_get(dim_p,dim_p);
	sf_point->pivot = px_get(dim_p);
	sf_point->LU_A = m_get(dim_p,dim_p);
	sf_point->gamma = v_get(dim_p);
	sf_point->inter = v_get(12);


	if (compute == 2)
	{

		sf_point->v_inter_1 = v_get(10);
		sf_point->v_inter_2 = v_get(10);

		sf_point->dphi_dk = v_get(10);
		sf_point->dp_dk = v_get(dim_p);


		sf_point->dA_dk = malloc(dim*sizeof(MAT*));
		sf_point->dB_dk = malloc(dim*sizeof(MAT*));

		// LU decompostion variables
		sf_point->gamma_k = v_get(dim_p);
		sf_point->RHS = v_get(dim_p);

		// Derivative variables
		for ( int k = 0 ; k < dim ; k++)
		{	
			sf_point->dA_dk[k] = m_get(dim_p,dim_p);
			sf_point->dB_dk[k] = m_get(dim_p,10);

		}


	}else if ( compute == 3)
	{

	}


	return sf_point;

}



shape_function * mls_shapefunction_materialpoint(MATERIAL_POINT * MP, int compute, MAT * nodes)
{


	// nodes and domain size from meshfree structure
	char * basis_type = "linear";
	char * weight = "cubic_spline";
	int dim = nodes->n;

	shape_function * sf_point = MP->shape_function;
	double * compute_point = MP->coords_n_1;
	int i ;



	// Initialise some storage matricies
	if ( sf_point == NULL){
		sf_point = new_shape_function(compute,basis_type,dim);
	}

	double xS[3] = {0,0,0};


	int dim_p = sf_point->B->m;

	switch(compute)
	{
		case(1):
		{
			// get point coordinates
			double * x = compute_point;


			// Find point neighbors ( have to do something smart with this)
			//sf_point->neighbours = get_materialpoint_neighbours(sf_point->neighbours,MP,nodes);
			int num_neighbours = MP->num_neighbours;


			// resize output vectors
			sf_point->phi = v_resize(sf_point->phi,num_neighbours);


			// get polynomial basis p(x), and derivatives if required
			polynomial_basis(sf_point->basis, x, dim, basis_type , compute);
			sf_point->p = get_col(sf_point->basis, 0, sf_point->p);


			// resize B
			if ( sf_point->B->n != num_neighbours)
			{
				sf_point->B = m_resize(sf_point->B, sf_point->B->m, num_neighbours);
			}

			// zero A and B
			m_zero(sf_point->B);
			sf_point->A = m_zero(sf_point->A);



			// construct moment matrix A, and B matrix
			for ( int j = 0 ; j < num_neighbours ; ++j)
			{

				// find domain size, and coords of node xi;
				double * xi = nodes->me[MP->neighbours->ive[j]];


				for ( int k = 0 ; k < dim ; k++)
				{

					xS[k] = x[k] - xi[k];
				}

				// get p(x_i)
				polynomial_basis(sf_point->basis_xi, xi, dim, basis_type, 1);
				sf_point->p_xi = get_col(sf_point->basis_xi,0,sf_point->p_xi);

				// get weight function for node I
				weight_function_materialpoint(sf_point->weights, MP, xS, dim);

				// shape function matricies
				// B;
				sf_point->bi = sv_mlt(sf_point->weights->ve[0],sf_point->p_xi, sf_point->bi);
				sf_point->B = set_col(sf_point->B,j,sf_point->bi);
				//A
				sf_point->ppT = v_outer_product(sf_point->p_xi, sf_point->p_xi, sf_point->ppT);
				sf_point->A =ms_mltadd(sf_point->A, sf_point->ppT, sf_point->weights->ve[0], sf_point->A);


			} // end loop over neighbours


			sf_point->LU_A = m_copy(sf_point->A,sf_point->LU_A);


			// lu decomposition of A
			tracecatch(
				LUfactor(sf_point->LU_A,sf_point->pivot);
				LUsolve(sf_point->LU_A,sf_point->pivot,sf_point->p,sf_point->gamma); 
				,"Shape_function");



			sf_point->phi = vm_mlt(sf_point->B,sf_point->gamma,sf_point->phi);

			double phi_sum = 0;
			phi_sum = v_sum(sf_point->phi);

			double tol = 1e-6;

			if ( fabs(phi_sum -1) > tol)
			{	

				printf("ERROR ERROR : PARTION OF UNITY FAILED \n");
				printf("phi_sum = %lf\n",phi_sum);
				v_foutput(stdout,sf_point->phi);

			}

			break;


} // end if compute == 1

case(2):
{



	// get point coordinates
	double * x = compute_point;

	// find neighbours of point x;
	//sf_point->neighbours = get_materialpoint_neighbours(sf_point->neighbours,MP,nodes);
	int num_neighbours = MP->num_neighbours;


	
	// get output sizes
	sf_point->phi = v_resize(sf_point->phi,num_neighbours);
	sf_point->dphi = m_resize(sf_point->dphi,num_neighbours,dim);


	// get polynomial basis p, and derivatives if required
	polynomial_basis(sf_point->basis, x, dim, basis_type , compute);
	sf_point->p = get_col(sf_point->basis, 0,sf_point->p);


	// resize B
	if ( sf_point->B->n != num_neighbours)
	{
		sf_point->B = m_resize(sf_point->B, sf_point->B->m, num_neighbours);
	}

	// zero A and B
	m_zero(sf_point->B);
	sf_point->A = m_zero(sf_point->A);


	// resize dB_dk and zero dA_dk
	for ( int k = 0 ; k < dim ; k++)
	{	
		if ( sf_point->dB_dk[k]->n != num_neighbours)
		{
			sf_point->dB_dk[k] = m_resize(sf_point->dB_dk[k], dim_p, num_neighbours);
			sf_point->dB_dk[k] = m_zero(sf_point->dB_dk[k]);

		}else{
			sf_point->dB_dk[k] = m_zero(sf_point->dB_dk[k]);
		}
		sf_point->dA_dk[k] = m_zero(sf_point->dA_dk[k]);
	}




	// construct moment matrix A, and B matrix
	for ( int j = 0 ; j < num_neighbours ; ++j)
	{
		// find domain size, and coords of node xi;
		double * xi = nodes->me[MP->neighbours->ive[j]];

		for ( int k = 0 ; k < dim ; k++)
		{
			xS[k] = x[k] - xi[k];
		}


		// get p(x_i)
		polynomial_basis(sf_point->basis_xi, xi, dim, basis_type , 1);
		sf_point->p_xi = get_col(sf_point->basis_xi,0,sf_point->p_xi);


		// get weight function for node I
		weight_function_materialpoint(sf_point->weights, MP, xS, dim);

		// shape function matricies
		// B;
		sf_point->bi = sv_mlt(sf_point->weights->ve[0], sf_point->p_xi, sf_point->bi);
		sf_point->B = set_col(sf_point->B,j,sf_point->bi);
				//A
		sf_point->ppT = v_outer_product(sf_point->p_xi, sf_point->p_xi, sf_point->ppT);
		sf_point->A = ms_mltadd(sf_point->A, sf_point->ppT, sf_point->weights->ve[0], sf_point->A);

		// shape function derivative matrix 
		for ( int k = 0 ; k < dim ; k++)
		{
			sf_point->dA_dk[k] = ms_mltadd(sf_point->dA_dk[k], sf_point->ppT, sf_point->weights->ve[1+k], sf_point->dA_dk[k]);
			sf_point->bi = sv_mlt(sf_point->weights->ve[1+k],sf_point->p_xi, sf_point->bi);
			sf_point->dB_dk[k] = set_col(sf_point->dB_dk[k],j,sf_point->bi);
		}



	} // end loop over neighbours


	// store intermediate results
	sf_point->v_inter_1 = v_resize(sf_point->v_inter_1, num_neighbours);
	sf_point->v_inter_2 = v_resize(sf_point->v_inter_2, num_neighbours);
	sf_point->LU_A = m_copy(sf_point->A,sf_point->LU_A);


	// lu decomposition of A
	tracecatch(
		LUfactor(sf_point->LU_A,sf_point->pivot);
		LUsolve(sf_point->LU_A,sf_point->pivot,sf_point->p,sf_point->gamma); 
		,"Shape_function");

	sf_point->phi = vm_mlt(sf_point->B,sf_point->gamma,sf_point->phi);
	
	double lc_sum[dim];
	// shape function derivatives
	for ( int k = 0; k < dim ; k++)
	{
		mv_mlt(sf_point->dA_dk[k],sf_point->gamma,sf_point->v_inter_1);
		get_col(sf_point->basis, k+1,sf_point->dp_dk);
		v_sub(sf_point->dp_dk,sf_point->v_inter_1,sf_point->RHS);
		LUsolve(sf_point->LU_A,sf_point->pivot,sf_point->RHS,sf_point->gamma_k);
		vm_mlt(sf_point->B,sf_point->gamma_k,sf_point->v_inter_1);
		vm_mlt(sf_point->dB_dk[k],sf_point->gamma,sf_point->v_inter_2);
		v_add(sf_point->v_inter_1, sf_point->v_inter_2, sf_point->dphi_dk);
		set_col(sf_point->dphi, k, sf_point->dphi_dk);
		lc_sum[k] = 0;
	}


	double phi_sum = 0;
	phi_sum = v_sum(sf_point->phi);
	double tol = 1e-4;

	assert(fabs(phi_sum -1) < tol);
	int index = 0;


	double maxphi = v_max(sf_point->phi,&index);

	if ( maxphi > 1)

	{	
		iv_foutput(stdout, MP->neighbours);
		assert(v_max(sf_point->phi,&index) < 1.00);

	}
	assert(v_max(sf_point->phi,&index) < 1.00);


	if ( fabs(phi_sum -1) > tol)
		{	

			printf("ERROR ERROR : PARTION OF UNITY FAILED \n");
			printf("ERROR ERROR : PARTION OF UNITY FAILED  \n");

		}


	// derivative consistencey check
	for ( int k = 0 ; k < num_neighbours; k++)
	{	
		double * xi = nodes->me[MP->neighbours->ive[k]];
		for ( int l = 0 ; l < dim ; l++)
		{
			lc_sum[l] += xi[l]*sf_point->dphi->me[k][l];
		}


	}
	// derivative consistencey check
	for ( int k = 0 ; k < dim ; k++)
		{
		if ( fabs(lc_sum[k]-1) > tol)
			{	
				printf("ERROR ERROR : FIRST ORDER DERIVATIVE CONSISTENCEY FAILED \n");
				printf("ERROR ERROR : FIRST ORDER DERIVATIVE CONSISTENCEY FAILED \n");

			}
		}






		break;




	} // end if compute == 2

	case(3):
	{
		printf("not yet supported");

		break;
	} // end loop over sample points


}// end of switch
	return sf_point;


}

