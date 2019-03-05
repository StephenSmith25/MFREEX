/* CREATE MASS MATRIX*/
#include "Integration/mass_matrix.h"

VEC *  mass_vector(MATERIAL_POINTS * MPS, meshfreeDomain * mfree){


	/* ------------------------------------------*/
	/* ----------------Mass Vector---------------*/
	/* ------------------------------------------*/

	VEC * phi;
	IVEC * neighbours;
	VEC * nodal_mass = v_get(mfree->num_nodes);

	int IS_AXI = MPS->IS_AXI;



	for ( int i = 0 ; i < MPS->num_material_points ; i++)
	{
		phi = MPS->sf_material_points->sf_list[i]->phi;
		neighbours = MPS->sf_material_points->sf_list[i]->neighbours;

		double volume = MPS->MP[i]->volume;
		double rho = MPS->MP[i]->rho;




		for ( int k = 0 ; k < neighbours->max_dim ; k++)
		{
			int index = neighbours->ive[k];
			nodal_mass->ve[index] += phi->ve[k]*rho*volume;
		}

	}



	return nodal_mass;
}