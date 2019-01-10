#include "Boundary/Displacement/setUpBC.h"





int setUpBC(EBC * ebc_, VEC * invMass, meshfreeDomain * mfree){
	


	// Initialsie matricies
	const int num_node_B = ebc_->nodes->max_dim;
	const int num_nodes = mfree->nodes->m;
	MAT * point = m_get(1,2);
	ebc_->V = m_get(num_nodes,num_node_B);
	ebc_->phi = m_get(num_node_B,num_nodes);
	MAT * invM = m_get(num_nodes,num_nodes);

	for ( int i = 0 ; i < num_nodes ; i++)
	{
		invM->me[i][i] = invMass->ve[i];
	}

	/*  Set coordinates */
	ebc_->coords = m_get(ebc_->nodes->max_dim,2);
	for ( int i = 0 ; i < ebc_->nodes->max_dim ;i++){
		ebc_->coords->me[i][0] = mfree->nodes->me[ebc_->nodes->ive[i]][0];
		ebc_->coords->me[i][1] = mfree->nodes->me[ebc_->nodes->ive[i]][1];
	}

	shape_function_container * sf_nodes = mls_shapefunction(ebc_->coords, "linear", "cubic", 2, 1, mfree);

	VEC * phi;
	IVEC * index;
	/*  Find v and phi, such that d_boundary =phi*{u} and Febc=V*Tebc */
	for ( int i = 0 ; i < ebc_->nodes->max_dim; i++){
		index = sf_nodes->sf_list[i]->neighbours;
		phi = sf_nodes->sf_list[i]->phi;

		for ( int m = 0 ; m < index->max_dim; m++){
			ebc_->V->me[index->ive[m]][i] += phi->ve[m] ; 
			ebc_->phi->me[i][index->ive[m]] += phi->ve[m];

		}
	}


	/*  Find P, inv(M)*V*(Phi*inv(M)*V)^-1 */
	ebc_->P= m_get(num_nodes,num_node_B);
	MAT * P_3_a = m_get(num_nodes,ebc_->nodes->max_dim);
	m_mlt(invM,ebc_->V,P_3_a);
	MAT * P_3_b = m_get(num_nodes,num_nodes);
	m_mlt(ebc_->phi,P_3_a,P_3_b);
	MAT * invP_3_b = m_get(ebc_->nodes->max_dim,ebc_->nodes->max_dim);
	m_inverse(P_3_b,invP_3_b); 
	MAT * P_3_c = m_get(num_nodes,ebc_->nodes->max_dim);
	m_mlt(ebc_->V,invP_3_b,P_3_c);
	m_mlt(invM,P_3_c,ebc_->P);	

	/*  Initialise  */
	ebc_->uCorrect1 = v_get(ebc_->nodes->max_dim);
	ebc_->uCorrect2 = v_get(ebc_->nodes->max_dim);

	// free any allocated memory that is not required anymore. 
	m_free(invM);
	m_free(point);
	m_free(P_3_a);
	m_free(P_3_b);
	m_free(P_3_c);
	m_free(invP_3_b);


	/*  Free the shape function storage phi */
	free_shapefunction_container(sf_nodes);


	return 0;
}
