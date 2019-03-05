
#include "mls_shapefunction.h"


shape_function_container *  mls_shapefunction(MAT * compute_points, int compute, meshfreeDomain * mfree)
{


	// nodes and domain size from meshfree structure



	int numnodes = mfree->num_nodes;
	MAT * nodes = mfree->nodes;
	char * shape = mfree->kernel_shape;
	char * basis_type = mfree->basis_type;
	char * weight = mfree->weight_function;
	int dim = mfree->dim;

	enum SUPPORT_TYPE support = mfree->kernel_support;


	VEC * domainSize = mfree->di;


	// number of sample points
	int numPoints = compute_points->m;

	shape_function ** sf_array = malloc( (numPoints) * sizeof(shape_function*) );
	shape_function_container * sf_list = malloc( 1 * sizeof(shape_function_container));


	for ( int i = 0 ; i < numPoints ; ++i)
	{
		sf_array[i] = malloc(1*sizeof(shape_function));
	}

	// test to ensure the program is running on a suitable number of threads
	if ( numPoints < 100)
	{
		omp_set_num_threads(1);
	}

	else if ( numPoints < 1000)
	{
		omp_set_num_threads(1); 
	}else if ( numPoints < 20000){

		omp_set_num_threads(4);

	}else{
		omp_set_num_threads(4);
	} 

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

	if ( compute == 1)
	{
		int i ;

#pragma omp parallel 
		{
			// initialise shape function calculation matricies (these will be free'd
			MAT * A = m_get(dim_p,dim_p);
			MAT * B = m_get(dim_p,12);
			MAT * basis_xi = m_get(dim_p,1);
			VEC * p_xi = v_get(dim_p);
			MAT * ppT = m_get(dim_p,dim_p);
			VEC * bi = v_get(dim_p);
			VEC * weights = v_get(1);
			MAT * invA = m_get(dim_p,dim_p);
			VEC * inter = v_get(12);
			MAT * basis = m_get(dim_p,1);
			VEC * p = v_get(dim_p);

			//printf("obtaining shape functions using %d threads\n", omp_get_num_threads());


			PERM * pivot = px_get(dim_p);
			MAT * LU_A = m_get(dim_p,dim_p);
			VEC * gamma = v_get(dim_p);


#pragma omp for 
			for (i = 0; i < numPoints; ++i){

			// get point coordinates
				double * x = compute_points->me[i];
			// initialise shape function structure
			// find neighbours of point x;
				IVEC * neighbours = point_neighbours(x,mfree);
				int num_neighbours = neighbours->max_dim;


			// get output sizes
				VEC * phi = v_get(num_neighbours);
				inter = v_resize(inter,num_neighbours);


			// get polynomial basis p, and derivatives if required
				polynomial_basis(basis, x, dim, basis_type , compute);
				p = get_col(basis, 0,p);


			// resize B
				if ( B->n != num_neighbours)
				{
					B = m_resize(B, dim_p, num_neighbours);
				}
				A = m_zero(A);


				// shifted coordinates, xi ( coords of point i) and domain size di
				double xS[3] = {0,0,0};

				// construct moment matrix A, and B matrix
				for ( int j = 0 ; j < num_neighbours ; ++j)
				{

					// find domain size, and coords of node xi;
					double * xi = nodes->me[neighbours->ive[j]];


					for ( int k = 0 ; k < dim ; k++)
					{

						xS[k] = x[k] - xi[k];
					}


					// get p(x_i)
					polynomial_basis(basis_xi, xi, dim, basis_type, 1);
					p_xi = get_col(basis_xi,0,p_xi);


					// get weight function for node I
					weight_function(weights, xS, neighbours->ive[j], mfree, compute);



					// shape function matricies
					// B;
					bi = sv_mlt(weights->ve[0], p_xi, bi);
					B = set_col(B,j,bi);
					//A
					ppT = v_outer_product(p_xi, p_xi, ppT);
					A =ms_mltadd(A, ppT, weights->ve[0], A);



			} // end loop over neighbours


			LU_A = m_copy(A,LU_A);


			// lu decomposition of A
			tracecatch(
				LUfactor(LU_A,pivot);
				LUsolve(LU_A,pivot,p,gamma); 
				,"Shape_function");



			// // invert A
			// invA = m_inverse(A,invA);
			// // get phi = p^T * inv(A) * B
			// inter = vm_mlt(invA, p, inter);
			// phi = vm_mlt(B,inter,phi);


			phi = vm_mlt(B,gamma,phi);

			double phi_sum = 0;
			phi_sum = v_sum(phi);

			double tol = 1e-6;

			if ( fabs(phi_sum -1) > tol)
			{	

				printf("ERROR ERROR : PARTION OF UNITY FAILED \n");
				printf("phi_sum = %lf\n",phi_sum);
				v_foutput(stdout,phi);


			}

			// set contents of shape function array 
			sf_array[i]->phi = phi;
			sf_array[i]->neighbours = neighbours;




		} // end loop over sample points

			// free memory
		M_FREE(A);
		M_FREE(B);
		M_FREE(basis_xi);
		V_FREE(p_xi);
		M_FREE(ppT);
		V_FREE(bi);
		V_FREE(weights);
		M_FREE(basis);
		V_FREE(p);
		M_FREE(invA);
		V_FREE(inter);


		PX_FREE(pivot);
		M_FREE(LU_A);
		V_FREE(gamma);


	} // end of loop over parallel region
} // end if compute == 1


if ( compute == 2)
{

#pragma omp parallel
	{


		// initialise shape function calculation matricies (these will be free'd)
		MAT * A = m_get(dim_p,dim_p);
		MAT * B = m_get(dim_p,10);
		MAT * basis_xi = m_get(dim_p,1);
		VEC * p_xi = v_get(dim_p);
		MAT * ppT = m_get(dim_p,dim_p);
		VEC * bi = v_get(dim_p);
		VEC * weights = v_get(1+dim);
		MAT * basis = m_get(1,1);
		VEC * p = v_get(dim_p);
		MAT * invA = m_get(dim_p,dim_p);


		// get derivatives of phi
		VEC * v_inter_1 = v_get(10);
		VEC * v_inter_2 = v_get(10);
		VEC * dphi = v_get(10);
		VEC * dp_dk = v_get(dim_p);


		MAT ** dA_dk = malloc(dim*sizeof(MAT*));
		MAT ** dB_dk = malloc(dim*sizeof(MAT*));


		// LU decompostion variables
		PERM * pivot = px_get(dim_p);
		MAT * LU_A = m_get(dim_p,dim_p);
		VEC * gamma = v_get(dim_p);
		VEC * gamma_k = v_get(dim_p);
		VEC * RHS = v_get(dim_p);

		// loop over to form derivative matricies
		for ( int k = 0 ; k < dim ; k++)
		{	
			dA_dk[k] = m_get(dim_p,dim_p);
			dB_dk[k] = m_get(dim_p,10);

		}

#pragma omp for 
		for (int i = 0; i < numPoints; ++i){


			// get point coordinates
			double * x = compute_points->me[i];

			// find neighbours of point x;
			IVEC * neighbours = point_neighbours(x,mfree);
			int num_neighbours = neighbours->max_dim;


			// get output sizes
			VEC * phi = v_get(num_neighbours);
			MAT * dphi_dk = m_get(num_neighbours,dim);


			// get polynomial basis p, and derivatives if required
			polynomial_basis(basis, x, dim, basis_type , compute);
			p = get_col(basis, 0,p);


			// resize B
			if ( B->n != num_neighbours)
			{
				B = m_resize(B, dim_p, num_neighbours);
			}
			A = m_zero(A);


			// resize dB_dk and zero dA_dk
			for ( int k = 0 ; k < dim ; k++)
			{	
				if ( dB_dk[k]->n != num_neighbours)
				{
					dB_dk[k] = m_resize(dB_dk[k], dim_p, num_neighbours);
					dB_dk[k] = m_zero(dB_dk[k]);

				}else{
					dB_dk[k] = m_zero(dB_dk[k]);
				}
				dA_dk[k] = m_zero(dA_dk[k]);
			}


			// shifted coordinates, xi ( coords of point i) and domain size di
			double xS[3] = {0,0,0};


			// construct moment matrix A, and B matrix
			for ( int j = 0 ; j < num_neighbours ; ++j)
			{

				// find domain size, and coords of node xi;
				double * xi = mfree->nodes->me[neighbours->ive[j]];

				for ( int k = 0 ; k < dim ; k++)
				{
						xS[k] = x[k] - xi[k];
				}

				// get p(x_i)
				polynomial_basis(basis_xi, xi, dim, basis_type , 1);
				p_xi = get_col(basis_xi,0,p_xi);


				// get weight function for node I
				weight_function(weights, xS, neighbours->ive[j], mfree, compute);

				// shape function matricies
				// B;
				bi = sv_mlt(weights->ve[0], p_xi, bi);
				B = set_col(B,j,bi);
				//A
				ppT = v_outer_product(p_xi, p_xi, ppT);
				A = ms_mltadd(A, ppT, weights->ve[0], A);

				// shape function derivative matrix 
				for ( int k = 0 ; k < dim ; k++)
				{
					dA_dk[k] = ms_mltadd(dA_dk[k], ppT, weights->ve[1+k], dA_dk[k]);
					bi = sv_mlt(weights->ve[1+k], p_xi, bi);
					dB_dk[k] = set_col(dB_dk[k],j,bi);
				}



			} // end loop over neighbours


			// store intermediate results
			v_inter_1 = v_resize(v_inter_1, num_neighbours);
			v_inter_2 = v_resize(v_inter_2, num_neighbours);

			LU_A = m_copy(A,LU_A);


			// lu decomposition of A
			tracecatch(
				LUfactor(LU_A,pivot);
				LUsolve(LU_A,pivot,p,gamma); 
				,"Shape_function");


			phi = vm_mlt(B,gamma,phi);

			double lc_sum[dim];
			// shape function derivatives
			for ( int k = 0; k < dim ; k++)
			{
				mv_mlt(dA_dk[k],gamma,v_inter_1);
				get_col(basis, k+1,dp_dk);
				v_sub(dp_dk,v_inter_1,RHS);
				LUsolve(LU_A,pivot,RHS,gamma_k);
				vm_mlt(B,gamma_k,v_inter_1);
				vm_mlt(dB_dk[k],gamma,v_inter_2);
				v_add(v_inter_1, v_inter_2, dphi);
				set_col(dphi_dk, k, dphi);
				lc_sum[k] = 0;
			}


			double phi_sum = 0;
			phi_sum = v_sum(phi);
			double tol = 1e-4;

			if ( fabs(phi_sum -1) > tol)
			{	

				printf("ERROR ERROR : PARTION OF UNITY FAILED  at node %d \n", i);
				printf("ERROR ERROR : PARTION OF UNITY FAILED  at node %d \n", i);

			}


			// derivative consistencey check
			for ( int k = 0 ; k < num_neighbours; k++)
			{	
				double * xi = mfree->nodes->me[neighbours->ive[k]];
				
				for ( int l = 0 ; l < dim ; l++)
				{
					lc_sum[l] += xi[l]*dphi_dk->me[k][l];
				}


			}
			// derivative consistencey check
			for ( int k = 0 ; k < dim ; k++)
			{
				if ( fabs(lc_sum[k]-1) > tol)
				{	
					printf("ERROR ERROR : FIRST ORDER DERIVATIVE CONSISTENCEY FAILED at node %d \n", i);
					printf("ERROR ERROR : FIRST ORDER DERIVATIVE CONSISTENCEY FAILED at node %d \n", i);
				}
			}


			// set contents of shape function array 
			sf_array[i]->phi = phi;
			sf_array[i]->dphi = dphi_dk;
			sf_array[i]->neighbours = neighbours;

			// assertion to catch any errors




		} // end loop over sample points

		// free memory
		// shape function matricies
		M_FREE(A);
		M_FREE(B);
		M_FREE(basis_xi);
		V_FREE(p_xi);
		M_FREE(ppT);
		V_FREE(bi);
		V_FREE(weights);
		M_FREE(basis);
		V_FREE(p);
		M_FREE(invA);


		// shape function derivative matricies
		V_FREE(v_inter_1);
		V_FREE(v_inter_2);
		V_FREE(dphi);
		V_FREE(dp_dk);

		// lu decompositon variables
		PX_FREE(pivot);
		M_FREE(LU_A);
		V_FREE(gamma);
		V_FREE(gamma_k);
		V_FREE(RHS);


		for ( int k = 0 ; k < dim ; k++)
		{
			M_FREE(dA_dk[k]);
			free(dA_dk[k]);
			M_FREE(dB_dk[k]);
		}
		free(dA_dk);
		free(dB_dk);





	}// end of omp parallel region

	} // end if compute == 2


	if ( compute == 3)
	{

#pragma omp parallel
		{


		// initialise shape function calculation matricies (these will be free'd)
			MAT * A = m_get(dim_p,dim_p);
			MAT * B = m_get(dim_p,10);
			MAT * basis_xi = m_get(dim_p,1);
			VEC * p_xi = v_get(dim_p);
			MAT * ppT = m_get(dim_p,dim_p);
			VEC * bi = v_get(dim_p);
			VEC * weights = v_get(1+dim);
			MAT * basis = m_get(1,1);
			VEC * p = v_get(dim_p);
			MAT * invA = m_get(dim_p,dim_p);


		// get derivatives of phi
			VEC * v_inter_1 = v_get(10);
			VEC * v_inter_2 = v_get(10);
			VEC * dphi = v_get(10);
			VEC * dp_dk = v_get(dim_p);
			MAT ** dA_dk = malloc(dim*sizeof(MAT*));
			MAT ** dB_dk = malloc(dim*sizeof(MAT*));

		// get second derivatives of phi
			MAT ** d2A_dk_dj = malloc(dim*dim*sizeof(MAT*));
			MAT ** d2B_dk_dj = malloc(dim*dim*sizeof(MAT*));
			VEC * d2phi = v_get(10);
			VEC * gamma_kj = v_get(dim_p);
			MAT * gamma_k_m = m_get(dim_p,dim);
			VEC * v_inter_3 = v_get(dim_p);

		// LU decompostion variables
			PERM * pivot = px_get(dim_p);
			MAT * LU_A = m_get(dim_p,dim_p);
			VEC * gamma = v_get(dim_p);
			VEC * gamma_k = v_get(dim_p);
			VEC * RHS = v_get(dim_p);



		// loop over to form derivative matricies
			for ( int k = 0 ; k < dim ; k++)
			{	
				dA_dk[k] = m_get(dim_p,dim_p);
				dB_dk[k] = m_get(dim_p,10);

				for ( int j = 0 ; j < dim ; ++j)
				{
					int indx = j + k*dim;
					d2A_dk_dj[indx] = m_get(dim_p,dim_p);
					d2B_dk_dj[indx] = m_get(dim_p,10);
				}

			}
#pragma omp for 
			for (int i = 0; i < numPoints; ++i){


			// get point coordinates
				double * x = compute_points->me[i];

			// initialise shape function structure
			// find neighbours of point x;
			IVEC * neighbours = point_neighbours(x,mfree);
				int num_neighbours = neighbours->max_dim;


			// get output sizes
				VEC * phi = v_get(num_neighbours);
				MAT * dphi_dk = m_get(num_neighbours,dim);
				MAT * d2phi_dk_dj = m_get(num_neighbours,dim*dim);


			// get polynomial basis p, and derivatives if required
				polynomial_basis(basis, x, dim, basis_type , compute);
				p = get_col(basis, 0,p);


			// resize B
				if ( B->n != num_neighbours)
				{
					B = m_resize(B, dim_p, num_neighbours);
				}
				A = m_zero(A);



			// resize dB_dk and zero dA_dk
				for ( int k = 0 ; k < dim ; k++)
				{	
					if ( dB_dk[k]->n != num_neighbours)
					{
						dB_dk[k] = m_resize(dB_dk[k], dim_p, num_neighbours);
						dB_dk[k] = m_zero(dB_dk[k]);

					}else{
						dB_dk[k] = m_zero(dB_dk[k]);
					}
					dA_dk[k] = m_zero(dA_dk[k]);

					for ( int j = 0 ; j < dim ; j++)
					{
						int indx = j + k*dim;
						if ( d2B_dk_dj[indx]->n != num_neighbours)
						{
							m_resize(d2B_dk_dj[indx],dim_p,num_neighbours);

						}
						m_zero(d2A_dk_dj[indx]);
					}
				}


			// shifted coordinates, xi ( coords of point i) and domain size di
				double xS[3] = {0,0,0};


			// construct moment matrix A, and B matrix
				for ( int j = 0 ; j < num_neighbours ; ++j)
				{

				// find domain size, and coords of node xi;
					double * xi = mfree->nodes->me[neighbours->ive[j]];
					

					for ( int k = 0 ; k < dim ; k++)
					{

						xS[k] = x[k] - xi[k];
					}


				// get p(x_i)
					polynomial_basis(basis_xi, xi, dim, basis_type , 1);
					p_xi = get_col(basis_xi,0,p_xi);


				// get weight function for node J
					weight_function(weights, xS, neighbours->ive[j], mfree, compute);



				// shape function matricies
				// B;
					bi = sv_mlt(weights->ve[0], p_xi, bi);
					B = set_col(B,j,bi);
				//A
					ppT = v_outer_product(p_xi, p_xi, ppT);
					A = ms_mltadd(A, ppT, weights->ve[0], A);

				// shape function derivative matrix 
					for ( int k = 0 ; k < dim ; k++)
					{
						dA_dk[k] = ms_mltadd(dA_dk[k], ppT, weights->ve[1+k], dA_dk[k]);
						bi = sv_mlt(weights->ve[1+k], p_xi, bi);
						dB_dk[k] = set_col(dB_dk[k],j,bi);

						for ( int l = 0 ; l < dim ; l++)
						{
							int indx = l + k*dim;
							d2A_dk_dj[indx] = ms_mltadd(d2A_dk_dj[indx],ppT,weights->ve[1+dim+indx],d2A_dk_dj[indx]);
							bi = sv_mlt(weights->ve[1+dim+indx],p_xi,bi);
							d2B_dk_dj[indx] = set_col(d2B_dk_dj[indx],j,bi);

						}
					}



			} // end loop over neighbours


			// store intermediate results
			v_inter_1 = v_resize(v_inter_1, num_neighbours);
			v_inter_2 = v_resize(v_inter_2, num_neighbours);

			LU_A = m_copy(A,LU_A);


			// lu decomposition of A
			tracecatch(
				LUfactor(LU_A,pivot);
				LUsolve(LU_A,pivot,p,gamma); 
				,"Shape_function");


			phi = vm_mlt(B,gamma,phi);

			double lc_sum[dim];
			// shape function derivatives
			for ( int k = 0; k < dim ; k++)
			{
				mv_mlt(dA_dk[k],gamma,v_inter_1);
				get_col(basis, k+1,dp_dk);
				v_sub(dp_dk,v_inter_1,RHS);
				LUsolve(LU_A,pivot,RHS,gamma_k);
				vm_mlt(B,gamma_k,v_inter_1);
				vm_mlt(dB_dk[k],gamma,v_inter_2);
				v_add(v_inter_1, v_inter_2, dphi);
				set_col(dphi_dk, k, dphi);


				// store gamma_k ( needed for shape function derivatives)
				set_col(gamma_k_m,k,gamma_k);


				lc_sum[k] = 0;


			}

			double qc_sum[dim*dim];


			for (int k = 0; k < dim; ++k)
			{


				for (int j = 0 ; j < dim ; j++)
				{
					int indx = j + (k*dim);
					// find RHS of lu decomposition such that LU gamma_j_k 
					
					// d2p_dk_dj
					get_col(basis,1+dim+indx,dp_dk);

					// dA_dk_dj * gamma
					mv_mlt(d2A_dk_dj[indx],gamma,v_inter_1);

					// dA_dj* gamma_k 
					gamma_k = get_col(gamma_k_m,k,gamma_k);
					mv_mlt(dA_dk[j],gamma_k,v_inter_2);

					// dA_dk * gamma_j
					gamma_k = get_col(gamma_k_m,j,gamma_k);
					mv_mlt(dA_dk[k],gamma_k,v_inter_3);

					// find RHS of the LU decomposition 
					v_sub(dp_dk,v_inter_1,RHS);
					v_sub(RHS,v_inter_2,RHS);
					v_sub(RHS,v_inter_3,RHS);

					LUsolve(LU_A,pivot,RHS,gamma_kj);

					// phi_k_j = gamma_k_j * B + gamma_j * B_k + gamma_k * B_j + gamma*B_k_j
					// gamma_k_j * B;
					vm_mlt(B,gamma_kj,v_inter_1);

					// gamma_j *B_k
					vm_mlt(dB_dk[k],gamma_k,v_inter_2);

					// gamma_k*B_j
					gamma_k = get_col(gamma_k_m,k,gamma_k);
					vm_mlt(dB_dk[j],gamma_k,v_inter_3);


					// gamma*B_k_j
					vm_mlt(d2B_dk_dj[indx],gamma,d2phi);

					v_add(d2phi,v_inter_1,d2phi);
					v_add(d2phi,v_inter_2,d2phi);
					v_add(d2phi,v_inter_3,d2phi);

					// set outputs
					set_col(d2phi_dk_dj, indx, d2phi);


					qc_sum[indx] = 0;




				}
			}
			double phi_sum = 0;
			phi_sum = v_sum(phi);
			double tol = 1e-4;

			if ( fabs(phi_sum -1) > tol)
			{	

				printf("ERROR ERROR : PARTION OF UNITY FAILED  at node %d \n", i);
				printf("ERROR ERROR : PARTION OF UNITY FAILED  at node %d \n", i);

			}


			// linear derivative consistencey check
			for ( int k = 0 ; k < num_neighbours; k++)
			{	
				double * xi = mfree->nodes->me[neighbours->ive[k]];
				
				for ( int l = 0 ; l < dim ; l++)
				{
					lc_sum[l] += xi[l]*dphi_dk->me[k][l];

					for ( int m = 0 ; m < dim ; m++)
					{
						int indx = m + l*dim;
						qc_sum[m + l*dim] += xi[l]*xi[m]*d2phi_dk_dj->me[k][indx];
					}
				}

			}

			// quadratic derivative consistencey check
			for ( int k = 0 ; k < dim ; k++)
			{
				if ( fabs(lc_sum[k]-1) > tol)
				{	
					printf("ERROR ERROR : FIRST ORDER DERIVATIVE CONSISTENCEY FAILED at node %d \n", i);
					printf("ERROR ERROR : FIRST ORDER DERIVATIVE CONSISTENCEY FAILED at node %d \n", i);
				}

				for ( int l = 0 ; l < dim ; l++)
				{
					int indx = l + k*dim;
					if ( l != k)
					{
						if ( fabs(qc_sum[indx] -1 ) > tol)
						{
							printf("ERROR ERROR : SECOND ORDER DERIVATIVE CONSISTENCEY FAILED at node %d \n", i);
							printf("ERROR ERROR : SECOND ORDER DERIVATIVE CONSISTENCEY FAILED at node %d \n", i);

						}
					}else{
						if ( fabs(qc_sum[indx] -2 ) > tol){
							printf("ERROR ERROR : SECOND ORDER DERIVATIVE CONSISTENCEY FAILED at node %d \n", i);
							printf("ERROR ERROR : SECOND ORDER DERIVATIVE CONSISTENCEY FAILED at node %d \n", i);

						}
					}
				}
			}


			// set contents of shape function array 
			sf_array[i]->phi = phi;
			sf_array[i]->dphi = dphi_dk;
			sf_array[i]->d2phi = d2phi_dk_dj;
			sf_array[i]->neighbours = neighbours;





		} // end loop over sample points

		// free memory
		// shape function matricies
		M_FREE(A);
		M_FREE(B);
		M_FREE(basis_xi);
		V_FREE(p_xi);
		M_FREE(ppT);
		V_FREE(bi);
		V_FREE(weights);
		M_FREE(basis);
		V_FREE(p);
		M_FREE(invA);


		// shape function derivative matricies
		V_FREE(v_inter_1);
		V_FREE(v_inter_2);
		V_FREE(dphi);
		V_FREE(dp_dk);


		// shape function second derivative matricies
		// get second derivatives of phi
		V_FREE(d2phi);
		V_FREE(gamma_kj);
		M_FREE(gamma_k_m);
		V_FREE(v_inter_3);

		// lu decompositon variables
		PX_FREE(pivot);
		M_FREE(LU_A);
		V_FREE(gamma);
		V_FREE(gamma_k);
		V_FREE(RHS);



		for ( int k = 0 ; k < dim ; k++)
		{
			M_FREE(dA_dk[k]);
			M_FREE(dB_dk[k]);

			for ( int j = 0 ; j < dim ; ++j)
			{
				int indx = j+ k*dim;
				M_FREE(d2A_dk_dj[indx]);
				M_FREE(d2B_dk_dj[indx]);
			}
		}
		free(dA_dk);
		free(dB_dk);
		free(d2A_dk_dj);
		free(d2B_dk_dj);




	}// end of omp parallel region

	} // end if compute == 3





	sf_list->number_of_points = numPoints;
	sf_list->compute = compute;
	sf_list->sf_list = sf_array;




	return sf_list;









}


int free_shapefunction_container(shape_function_container * sf_container)
{	

	for (int i = 0 ; i < sf_container->number_of_points ; i++)
	{
		V_FREE(sf_container->sf_list[i]->phi);
		IV_FREE(sf_container->sf_list[i]->neighbours);

		if ( sf_container->compute == 2)
		{
			M_FREE(sf_container->sf_list[i]->dphi);

		}

		if ( sf_container->compute == 3)
		{
			M_FREE(sf_container->sf_list[i]->dphi);
			M_FREE(sf_container->sf_list[i]->d2phi);

		}
		free(sf_container->sf_list[i]);
	}

	
	free(sf_container->sf_list);
	free(sf_container);
	return 0;
}



int update_mls_shapefunction(shape_function * sf_point, double compute_point[3], int compute, meshfreeDomain * mfree)
{


	// nodes and domain size from meshfree structure
	MAT * nodes = mfree->nodes;
	char * shape = mfree->kernel_shape;
	char * basis_type = mfree->basis_type;
	char * weight = mfree->weight_function;
	int dim = mfree->dim;
	enum SUPPORT_TYPE support = mfree->kernel_support;
	VEC * domainSize = mfree->di;



	int dim_p = sf_point->basis_xi->max_size;


	// initialise shape function calculation matricies (these will be free'd
	MAT * A = sf_point->A ;
	MAT * B = sf_point->B ;
	MAT * basis_xi = sf_point->basis_xi ;
	VEC * p_xi = sf_point->p_xi ;
	MAT * ppT = sf_point->ppT ;
	VEC * bi = sf_point->bi ;
	VEC * weights = sf_point->weights ;
	MAT * invA = sf_point->invA ;
	VEC * inter = sf_point->inter ;
	MAT * basis = sf_point->basis ;
	VEC * p = sf_point->p ;
	PERM * pivot = sf_point->pivot ;
	MAT * LU_A = sf_point->LU_A ;
	VEC * gamma = sf_point->gamma ;
	int i ;


	VEC * phi = sf_point->phi;
	MAT * dphi = sf_point->dphi;
	MAT * d2phi = sf_point->d2phi;


	iv_free(sf_point->neighbours);





	// Just first oder shape function
	if ( compute == 1)
	{
			// get point coordinates
			double * x = compute_point;



			// Find point neighbors ( have to do something smart with this)
			IVEC * neighbours = point_neighbours(x,mfree);
			int num_neighbours = neighbours->max_dim;


			// get output sizes
			phi = v_resize(phi,num_neighbours);
			inter = v_resize(inter,num_neighbours);


			// get polynomial basis p, and derivatives if required
			polynomial_basis(basis, x, dim, basis_type , compute);
			p = get_col(basis, 0,p);


			// resize B
			if ( B->n != num_neighbours)
			{
				B = m_resize(B, dim_p, num_neighbours);
			}
				A = m_zero(A);


			// shifted coordinates, xi ( coords of point i) and domain size di
			double xS[3] = {0,0,0};

			// construct moment matrix A, and B matrix
			for ( int j = 0 ; j < num_neighbours ; ++j)
				{

				// find domain size, and coords of node xi;
				double * xi = nodes->me[neighbours->ive[j]];
				for ( int k = 0 ; k < dim ; k++)
				{

					xS[k] = x[k] - xi[k];
				}

				// get p(x_i)
				polynomial_basis(basis_xi, xi, dim, basis_type, 1);
				p_xi = get_col(basis_xi,0,p_xi);

				// get weight function for node I
				weight_function(weights, xS, neighbours->ive[j], mfree, compute);

				// shape function matricies
				// B;
				bi = sv_mlt(weights->ve[0], p_xi, bi);
				B = set_col(B,j,bi);
				//A
				ppT = v_outer_product(p_xi, p_xi, ppT);
				A =ms_mltadd(A, ppT, weights->ve[0], A);


			} // end loop over neighbours


			LU_A = m_copy(A,LU_A);


			// lu decomposition of A
			tracecatch(
				LUfactor(LU_A,pivot);
				LUsolve(LU_A,pivot,p,gamma); 
				,"Shape_function");



			phi = vm_mlt(B,gamma,phi);

			double phi_sum = 0;
			phi_sum = v_sum(phi);

			double tol = 1e-6;

			if ( fabs(phi_sum -1) > tol)
			{	

				printf("ERROR ERROR : PARTION OF UNITY FAILED \n");
				printf("phi_sum = %lf\n",phi_sum);
				v_foutput(stdout,phi);

			}

			// set contents of shape function array 
			sf_point->phi = phi;
			sf_point->neighbours = neighbours;


} // end if compute == 1


if ( compute == 2)
{


			VEC * v_inter_1 = sf_point->v_inter_1 ;
			VEC * v_inter_2 = sf_point->v_inter_2 ;
			VEC * dphi_k = sf_point->dphia ;
			VEC * dp_dk = sf_point->dp_dk ;
			MAT ** dA_dk = sf_point->dA_dk ;
			MAT ** dB_dk = sf_point->dB_dk ;
			VEC * gamma_k = sf_point->gamma_k ;
			VEC * RHS = sf_point->RHS ;
			VEC * gamma = sf_point->gamma ;

			// get point coordinates
			double * x = compute_point;

			// find neighbours of point x;
			IVEC * neighbours = point_neighbours(x,mfree);
			int num_neighbours = neighbours->max_dim;


			// get output sizes
			phi = v_resize(phi,num_neighbours);
			dphi = m_resize(dphi,num_neighbours,dim);


			// get polynomial basis p, and derivatives if required
			polynomial_basis(basis, x, dim, basis_type , compute);
			p = get_col(basis, 0,p);


			// resize B
			if ( B->n != num_neighbours)
			{
				B = m_resize(B, dim_p, num_neighbours);
			}
			A = m_zero(A);


			// resize dB_dk and zero dA_dk
			for ( int k = 0 ; k < dim ; k++)
			{	
				if ( dB_dk[k]->n != num_neighbours)
				{
					dB_dk[k] = m_resize(dB_dk[k], dim_p, num_neighbours);
					dB_dk[k] = m_zero(dB_dk[k]);

				}else{
					dB_dk[k] = m_zero(dB_dk[k]);
				}
				dA_dk[k] = m_zero(dA_dk[k]);
			}


			// shifted coordinates, xi ( coords of point i) and domain size di
			double xS[3] = {0,0,0};


			// construct moment matrix A, and B matrix
			for ( int j = 0 ; j < num_neighbours ; ++j)
			{

				// find domain size, and coords of node xi;
				double * xi = mfree->nodes->me[neighbours->ive[j]];

				for ( int k = 0 ; k < dim ; k++)
				{
						xS[k] = x[k] - xi[k];
				}

				// get p(x_i)
				polynomial_basis(basis_xi, xi, dim, basis_type , 1);
				p_xi = get_col(basis_xi,0,p_xi);


				// get weight function for node I
				weight_function(weights, xS, neighbours->ive[j], mfree, compute);

				// shape function matricies
				// B;
				bi = sv_mlt(weights->ve[0], p_xi, bi);
				B = set_col(B,j,bi);
				//A
				ppT = v_outer_product(p_xi, p_xi, ppT);
				A = ms_mltadd(A, ppT, weights->ve[0], A);

				// shape function derivative matrix 
				for ( int k = 0 ; k < dim ; k++)
				{
					dA_dk[k] = ms_mltadd(dA_dk[k], ppT, weights->ve[1+k], dA_dk[k]);
					bi = sv_mlt(weights->ve[1+k], p_xi, bi);
					dB_dk[k] = set_col(dB_dk[k],j,bi);
				}



			} // end loop over neighbours


			// store intermediate results
			v_inter_1 = v_resize(v_inter_1, num_neighbours);
			v_inter_2 = v_resize(v_inter_2, num_neighbours);

			LU_A = m_copy(A,LU_A);


			// lu decomposition of A
			tracecatch(
				LUfactor(LU_A,pivot);
				LUsolve(LU_A,pivot,p,gamma); 
				,"Shape_function");


			phi = vm_mlt(B,gamma,phi);

			double lc_sum[dim];
			// shape function derivatives
			for ( int k = 0; k < dim ; k++)
			{
				mv_mlt(dA_dk[k],gamma,v_inter_1);
				get_col(basis, k+1,dp_dk);
				v_sub(dp_dk,v_inter_1,RHS);
				LUsolve(LU_A,pivot,RHS,gamma_k);
				vm_mlt(B,gamma_k,v_inter_1);
				vm_mlt(dB_dk[k],gamma,v_inter_2);
				v_add(v_inter_1, v_inter_2, dphi_k);
				set_col(dphi, k, dphi_k);
				lc_sum[k] = 0;
			}


			double phi_sum = 0;
			phi_sum = v_sum(phi);
			double tol = 1e-4;

			if ( fabs(phi_sum -1) > tol)
			{	

				printf("ERROR ERROR : PARTION OF UNITY FAILED \n");
				printf("ERROR ERROR : PARTION OF UNITY FAILED  \n");

			}


			// derivative consistencey check
			for ( int k = 0 ; k < num_neighbours; k++)
			{	
				double * xi = mfree->nodes->me[neighbours->ive[k]];
				
				for ( int l = 0 ; l < dim ; l++)
				{
					lc_sum[l] += xi[l]*dphi->me[k][l];
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

			sf_point->neighbours = neighbours;



	} // end if compute == 2


	if ( compute == 3)
	{
			// LU decompostion variables
			VEC * gamma_k = sf_point->gamma_k ;
			VEC * RHS = sf_point->RHS ;


			VEC * gamma = sf_point->gamma ;

			// get derivatives of phi
			VEC * v_inter_1 = sf_point->v_inter_1 ;
			VEC * v_inter_2 = sf_point->v_inter_2 ;
			VEC * dphi_dk = sf_point->dphia ;
			VEC * dp_dk = sf_point->dp_dk ;
			MAT ** dA_dk = sf_point->dA_dk ;
			MAT ** dB_dk = sf_point->dB_dk ;

			MAT ** d2A_dk_dj = sf_point->d2A_dk_dj ;
			MAT ** d2B_dk_dj = sf_point->d2B_dk_dj ;
			VEC * d2phi_dk_dj = sf_point->d2phia ;
			VEC * gamma_kj = sf_point->gamma_kj ;
			MAT * gamma_k_m = sf_point->gamma_k_m ;
			VEC * v_inter_3 = sf_point->v_inter_3 ;



			// get point coordinates
			double * x = compute_point;

			// initialise shape function structure
			// find neighbours of point x;
			IVEC * neighbours = point_neighbours(x,mfree);
			int num_neighbours = neighbours->max_dim;


			// get output sizes
			phi = v_resize(phi,num_neighbours);
			dphi = m_resize(dphi,num_neighbours,dim);
			d2phi = m_resize(d2phi,num_neighbours,dim*dim);


			// get polynomial basis p, and derivatives if required
			polynomial_basis(basis, x, dim, basis_type , compute);
			p = get_col(basis, 0,p);


			// resize B
				if ( B->n != num_neighbours)
				{
					B = m_resize(B, dim_p, num_neighbours);
				}
				A = m_zero(A);



			// resize dB_dk and zero dA_dk
				for ( int k = 0 ; k < dim ; k++)
				{	
					if ( dB_dk[k]->n != num_neighbours)
					{
						dB_dk[k] = m_resize(dB_dk[k], dim_p, num_neighbours);
						dB_dk[k] = m_zero(dB_dk[k]);

					}else{
						dB_dk[k] = m_zero(dB_dk[k]);
					}
					dA_dk[k] = m_zero(dA_dk[k]);

					for ( int j = 0 ; j < dim ; j++)
					{
						int indx = j + k*dim;
						if ( d2B_dk_dj[indx]->n != num_neighbours)
						{
							m_resize(d2B_dk_dj[indx],dim_p,num_neighbours);

						}
						m_zero(d2A_dk_dj[indx]);
					}
				}


			// shifted coordinates, xi ( coords of point i) and domain size di
				double xS[3] = {0,0,0};


			// construct moment matrix A, and B matrix
				for ( int j = 0 ; j < num_neighbours ; ++j)
				{

				// find domain size, and coords of node xi;
					double * xi = mfree->nodes->me[neighbours->ive[j]];
					

					for ( int k = 0 ; k < dim ; k++)
					{

						xS[k] = x[k] - xi[k];
					}


				// get p(x_i)
					polynomial_basis(basis_xi, xi, dim, basis_type , 1);
					p_xi = get_col(basis_xi,0,p_xi);


				// get weight function for node J
					weight_function(weights, xS, neighbours->ive[j], mfree, compute);



				// shape function matricies
				// B;
					bi = sv_mlt(weights->ve[0], p_xi, bi);
					B = set_col(B,j,bi);
				//A
					ppT = v_outer_product(p_xi, p_xi, ppT);
					A = ms_mltadd(A, ppT, weights->ve[0], A);

				// shape function derivative matrix 
					for ( int k = 0 ; k < dim ; k++)
					{
						dA_dk[k] = ms_mltadd(dA_dk[k], ppT, weights->ve[1+k], dA_dk[k]);
						bi = sv_mlt(weights->ve[1+k], p_xi, bi);
						dB_dk[k] = set_col(dB_dk[k],j,bi);

						for ( int l = 0 ; l < dim ; l++)
						{
							int indx = l + k*dim;
							d2A_dk_dj[indx] = ms_mltadd(d2A_dk_dj[indx],ppT,weights->ve[1+dim+indx],d2A_dk_dj[indx]);
							bi = sv_mlt(weights->ve[1+dim+indx],p_xi,bi);
							d2B_dk_dj[indx] = set_col(d2B_dk_dj[indx],j,bi);

						}
					}



			} // end loop over neighbours


			// store intermediate results
			v_inter_1 = v_resize(v_inter_1, num_neighbours);
			v_inter_2 = v_resize(v_inter_2, num_neighbours);

			LU_A = m_copy(A,LU_A);


			// lu decomposition of A
			tracecatch(
				LUfactor(LU_A,pivot);
				LUsolve(LU_A,pivot,p,gamma); 
				,"Shape_function");


			phi = vm_mlt(B,gamma,phi);

			double lc_sum[dim];
			// shape function derivatives
			for ( int k = 0; k < dim ; k++)
			{
				mv_mlt(dA_dk[k],gamma,v_inter_1);
				get_col(basis, k+1,dp_dk);
				v_sub(dp_dk,v_inter_1,RHS);
				LUsolve(LU_A,pivot,RHS,gamma_k);
				vm_mlt(B,gamma_k,v_inter_1);
				vm_mlt(dB_dk[k],gamma,v_inter_2);
				v_add(v_inter_1, v_inter_2, dphi_dk);
				set_col(dphi, k, dphi_dk);


				// store gamma_k ( needed for shape function derivatives)
				set_col(gamma_k_m,k,gamma_k);


				lc_sum[k] = 0;


			}

			double qc_sum[dim*dim];


			for (int k = 0; k < dim; ++k)
			{


				for (int j = 0 ; j < dim ; j++)
				{
					int indx = j + (k*dim);
					// find RHS of lu decomposition such that LU gamma_j_k 
					
					// d2p_dk_dj
					get_col(basis,1+dim+indx,dp_dk);

					// dA_dk_dj * gamma
					mv_mlt(d2A_dk_dj[indx],gamma,v_inter_1);

					// dA_dj* gamma_k 
					gamma_k = get_col(gamma_k_m,k,gamma_k);
					mv_mlt(dA_dk[j],gamma_k,v_inter_2);

					// dA_dk * gamma_j
					gamma_k = get_col(gamma_k_m,j,gamma_k);
					mv_mlt(dA_dk[k],gamma_k,v_inter_3);

					// find RHS of the LU decomposition 
					v_sub(dp_dk,v_inter_1,RHS);
					v_sub(RHS,v_inter_2,RHS);
					v_sub(RHS,v_inter_3,RHS);

					LUsolve(LU_A,pivot,RHS,gamma_kj);

					// phi_k_j = gamma_k_j * B + gamma_j * B_k + gamma_k * B_j + gamma*B_k_j
					// gamma_k_j * B;
					vm_mlt(B,gamma_kj,v_inter_1);

					// gamma_j *B_k
					vm_mlt(dB_dk[k],gamma_k,v_inter_2);

					// gamma_k*B_j
					gamma_k = get_col(gamma_k_m,k,gamma_k);
					vm_mlt(dB_dk[j],gamma_k,v_inter_3);


					// gamma*B_k_j
					vm_mlt(d2B_dk_dj[indx],gamma,d2phi_dk_dj);

					v_add(d2phi_dk_dj,v_inter_1,d2phi_dk_dj);
					v_add(d2phi_dk_dj,v_inter_2,d2phi_dk_dj);
					v_add(d2phi_dk_dj,v_inter_3,d2phi_dk_dj);

					// set outputs
					set_col(d2phi, indx, d2phi_dk_dj);


					qc_sum[indx] = 0;




				}
			}
			double phi_sum = 0;
			phi_sum = v_sum(phi);
			double tol = 1e-4;

			if ( fabs(phi_sum -1) > tol)
			{	

				printf("ERROR ERROR : PARTION OF UNITY FAILED  \n");
				printf("ERROR ERROR : PARTION OF UNITY FAILED \n");

			}


			// linear derivative consistencey check
			for ( int k = 0 ; k < num_neighbours; k++)
			{	
				double * xi = mfree->nodes->me[neighbours->ive[k]];
				
				for ( int l = 0 ; l < dim ; l++)
				{
					lc_sum[l] += xi[l]*dphi->me[k][l];

					for ( int m = 0 ; m < dim ; m++)
					{
						int indx = m + l*dim;
						qc_sum[m + l*dim] += xi[l]*xi[m]*d2phi->me[k][indx];
					}
				}

			}

			// quadratic derivative consistencey check
			for ( int k = 0 ; k < dim ; k++)
			{
				if ( fabs(lc_sum[k]-1) > tol)
				{	
					printf("ERROR ERROR : FIRST ORDER DERIVATIVE CONSISTENCEY FAILED  \n");
					printf("ERROR ERROR : FIRST ORDER DERIVATIVE CONSISTENCEY FAILED  \n");
				}

				for ( int l = 0 ; l < dim ; l++)
				{
					int indx = l + k*dim;
					if ( l != k)
					{
						if ( fabs(qc_sum[indx] -1 ) > tol)
						{
							printf("ERROR ERROR : SECOND ORDER DERIVATIVE CONSISTENCEY FAILED \n");
							printf("ERROR ERROR : SECOND ORDER DERIVATIVE CONSISTENCEY FAILED n");

						}
					}else{
						if ( fabs(qc_sum[indx] -2 ) > tol){
							printf("ERROR ERROR : SECOND ORDER DERIVATIVE CONSISTENCEY FAILED  \n");
							printf("ERROR ERROR : SECOND ORDER DERIVATIVE CONSISTENCEY FAILED  \n");

						}
					}
				}
			}


			// set contents of shape function array 
			sf_point->neighbours = neighbours;


	} // end loop over sample points



	return 0;

}


