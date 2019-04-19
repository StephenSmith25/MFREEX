#include "BLOCKS.h"

// Domain constructors  
int AddBlockSetToDomain(DOMAIN * domain, int ID )
{
	// create a new blockset with ID
	BLOCKSET * newBlockset = malloc(1*sizeof(BLOCKSET)) ;
	newBlockset->ID = ID;
	newBlockset->elements = NULL;
	// add to domain
	newBlockset->next = domain->blocksets;
	domain->blocksets = newBlockset;

	domain->NUM_BLOCK_SETS++;

	return 0;
}
int AddSideSetToDomain(DOMAIN * domain, int ID )
{
	// create a new sideset with ID
	SIDESET * newSideSet = malloc(1*sizeof(SIDESET)) ;
	newSideSet->ID = ID;
	// add to domain
	newSideSet->next = domain->sidesets;
	domain->sidesets = newSideSet;


	domain->NUM_SIDE_SETS++;

	return 0;
}
int AddNodeSetToDomain(DOMAIN * domain, int ID )
{
	// create a new nodeset with ID
	NODESET * newNodeSet = malloc(1*sizeof(NODESET));
	newNodeSet->ID = ID;
	newNodeSet->dof_constrained = NULL;
	// add to domain
	newNodeSet->next = domain->nodesets;
	domain->nodesets = newNodeSet;


	domain->NUM_NODE_SETS++;



	return 0;
}

int AddPhysicalNameToDomain(DOMAIN * domain, PHYSICAL_TYPE type, int ID)
{
	// create a new physical name with ID and Type
	PHYSICAL_NAME * physical_name = malloc(1*sizeof(PHYSICAL_NAME));
	physical_name->type = type;
	physical_name->ID = ID;
	
	// add to domain
	physical_name->next = domain->physical_names;
	domain->physical_names = physical_name;



	return 0;


}


PHYSICAL_NAME * FindPhysicalNameByID(DOMAIN * domain, int ID)
{
	PHYSICAL_NAME * current = domain->physical_names;
	int testID;
	while ( (current != NULL) && (testID != ID) )
	{
		testID = current->ID;
		if ( testID == ID)
			return current;
		// else move on
		current = current->next;
	}



	return NULL;
}

BLOCKSET * FindBlockSetByID(DOMAIN * domain, int ID)
{
	BLOCKSET * current = domain->blocksets;
	int testID = -1;

	while (( current != NULL) && ( testID != ID))
	{
		testID = current->ID;

		if ( testID == ID)
		{
			return current;
		}
		// else move on

		current = current->next;
	}

	printf("COULDNT FIND THE BLOCKSET WITH ID GIVEN \n");
	return NULL;



}

SIDESET * FindSideSetByID(DOMAIN * domain, int ID)
{


	SIDESET* current = domain->sidesets;
	int testID = -1;

	while (( current != NULL) && ( testID != ID))
	{
		testID = current->ID;

		if ( testID == ID)
		{
			return current;
		}
		// else move on

		current = current->next;
	}

	printf("COULDNT FIND THE SIDESET WITH ID GIVEN \n");
	return NULL;



}

NODESET * FindNodeSetByID(DOMAIN * domain, int ID)
{


	NODESET * current = domain->nodesets;
	int testID = -1;

	while (( current != NULL) && ( testID != ID))
	{
		testID = current->ID;

		if ( testID == ID)
		{
			return current;
		}
		// else move on

		current = current->next;
	}

	printf("COULDNT FIND THE NODESET WITH ID GIVEN \n");
	return NULL;



}



ELEMENT * AddElementToBlockSet(BLOCKSET * blockset, ELEMENT_TYPE etype, int * verticies)
{

	// Create a new element

	ELEMENT * newElement = CreateNewElement(etype, verticies);
	newElement->next = NULL;


	// Save refrence to the header to return
	ELEMENT * ret = blockset->elements;
	

	// add new element to tail
	ELEMENT * end = blockset->elements;

	if ( end == NULL)
	{
		return newElement;

	}else{

		// loop until find the end of list
		while ( end->next != NULL)
		{
			end = end->next;
		}
		end->next = newElement;
	}

	return ret;
}

ELEMENT * AddElementToSideSet(SIDESET * sideset, ELEMENT_TYPE etype, int * verticies)
{

	// Create a new element
	ELEMENT * newElement = CreateNewElement(etype, verticies);
	newElement->next = NULL;


	// Save refrence to the header to return
	ELEMENT * ret = sideset->elements;
	

	// add new element to tail
	ELEMENT * end = sideset->elements;

	if ( end == NULL)
	{
		return newElement;

	}else{

		// loop until find the end of list
		while ( end->next != NULL)
		{
			end = end->next;
		}
		end->next = newElement;
	}

	return ret;
}

ELEMENT * AddElementToNodeSet(NODESET * nodeset, ELEMENT_TYPE etype, int * verticies)
{

	// Create a new element
	ELEMENT * newElement = CreateNewElement(etype, verticies);
	newElement->next = NULL;


	// Save refrence to the header to return
	ELEMENT * ret = nodeset->elements;
	

	// add new element to tail
	ELEMENT * end = nodeset->elements;

	if ( end == NULL)
	{
		return newElement;

	}else{

		// loop until find the end of list
		while ( end->next != NULL)
		{
			end = end->next;
		}
		end->next = newElement;
	}

	return ret;
}


int WriteBlockSetElementsToFile(BLOCKSET * blockset, char * FILENAME)
{

	FILE * fp = fopen(FILENAME,"w");
	if (fp == NULL) {
   	 	perror("Failed: ");
    	return 1;
	}


	ELEMENT * current = blockset->elements;


	while (current != NULL)
	{
		if (current->etype == TRIANGLES)
			fprintf(fp,"%d,%d,%d\n",current->verticies[0],current->verticies[1],current->verticies[2]);
		else if ( current->etype == LINE)
			fprintf(fp,"%d,%d\n",current->verticies[0],current->verticies[1]);
		else if ( current->etype == QUAD)
			fprintf(fp,"%d,%d,%d,%d\n",current->verticies[0],current->verticies[1],
				current->verticies[2],current->verticies[3]);
		else if ( current->etype == TETRAHEDRON)
			fprintf(fp,"%d,%d,%d,%d\n",current->verticies[0],current->verticies[1],
				current->verticies[2],current->verticies[3]);
		else
			fprintf(stderr,"ELEMENT TYPE NOT SUPPORTED \n");

		// move to next element;
		current = current->next;
	}

	return 0;

}

int WriteSideSetElementsToFile(SIDESET * sideset, char * FILENAME)
{

	FILE * fp = fopen(FILENAME,"w");
	if (fp == NULL) {
   	 	perror("Failed: ");
    	return 1;
	}


	ELEMENT * current = sideset->elements;


	while (current != NULL)
	{
		if (current->etype == TRIANGLES)
			fprintf(fp,"%d,%d,%d\n",current->verticies[0],current->verticies[1],current->verticies[2]);
		else if ( current->etype == LINE)
			fprintf(fp,"%d,%d\n",current->verticies[0],current->verticies[1]);
		else if ( current->etype == QUAD)
			fprintf(fp,"%d,%d,%d,%d\n",current->verticies[0],current->verticies[1],
				current->verticies[2],current->verticies[3]);
		else if ( current->etype == TETRAHEDRON)
			fprintf(fp,"%d,%d,%d,%d\n",current->verticies[0],current->verticies[1],
				current->verticies[2],current->verticies[3]);
		else
			fprintf(stderr,"ELEMENT TYPE NOT SUPPORTED \n");

		// move to next element;
		current = current->next;
	}

	return 0;

}


int WriteNodeSetElementsToFile(NODESET * nodeset, char * FILENAME)
{

	FILE * fp = fopen(FILENAME,"w");
	if (fp == NULL) {
   	 	perror("Failed: ");
    	return 1;
	}


	ELEMENT * current = nodeset->elements;


	while (current != NULL)
	{
		if (current->etype == TRIANGLES)
			fprintf(fp,"%d,%d,%d\n",current->verticies[0],current->verticies[1],current->verticies[2]);
		else if ( current->etype == LINE)
			fprintf(fp,"%d,%d\n",current->verticies[0],current->verticies[1]);
		else if ( current->etype == QUAD)
			fprintf(fp,"%d,%d,%d,%d\n",current->verticies[0],current->verticies[1],
				current->verticies[2],current->verticies[3]);
		else if ( current->etype == TETRAHEDRON)
			fprintf(fp,"%d,%d,%d,%d\n",current->verticies[0],current->verticies[1],
				current->verticies[2],current->verticies[3]);
		else
			fprintf(stderr,"ELEMENT TYPE NOT SUPPORTED \n");

		// move to next element;
		current = current->next;
	}

	return 0;

}



//DisplacmentBCS;
int AddDOFConstraintToNodeSet(NODESET * nodeset,DOF_CONSTRAINT * new_constraint )
{

	// add new element to tail
	DOF_CONSTRAINT * end = nodeset->dof_constrained;

	if ( end == NULL)
	{
		nodeset->dof_constrained = new_constraint;
		return 0;

	}else{

		// loop until find the end of list
		while ( end->next != NULL)
		{
			end = end->next;
		}
		end->next = new_constraint;
	}


	return 0;


}

DOF_CONSTRAINT * FindDOFConstraintInNodeSet(NODESET * nodeset, DIRECTION dir)
{
	// loop over each constraint in nodeset to find one that adheres to the direction:
	DOF_CONSTRAINT * dof_constraint = nodeset->dof_constrained; 
	while (dof_constraint != NULL)
	{
		if (dof_constraint->dir ==dir)
		{
			return dof_constraint;
		}

		dof_constraint = dof_constraint->next;
	}


	return NULL;
}


#include "Integration/quad.h"
int GetMaterialPointsBlockSet(DOMAIN * domain , int order)
{

	int number_of_quadrature_points = 0;

	BLOCKSET * blockset = domain->blocksets;


	while ( blockset != NULL){


	//  create quadrature points
    // LOOP OVER EACH ELEMENT 
	ELEMENT * element;
	while ( element != NULL)
	{	
		// GET QUADRATURE POINTS AT EACH ELEMENT 
		QUAD * quad = GetQuadPoints(domain->NODES, element, order);
		element = element->next;

		// Create a material point for each quadarture-point
		for ( int i = 0 ; i < quad->points->m ; i++)
		{

		}


		// Add material point to blockset

		// prototype
		// AddMaterialPointToBlockSet(BLOCKSET * blocket, MATERIAL_POINT * mp)








	}

	blockset = blockset->next;
	}

	// 	QUAD_TRIANGLE * quad_triangles = create_triangle_quadrature_points(tri_points, 
	// 	triangles,number_of_triangles, quad_orders );



	// 	int number_of_quadrature_points = quad_triangles->QUAD_POINTS->m;
	// 	material_points->num_material_points = number_of_quadrature_points;


	// 	MATERIAL_POINT ** MPS = malloc(number_of_quadrature_points * sizeof(MATERIAL_POINT));



	// 	double point_coords[3] = {0,0,0};




	// 	for ( int i = 0 ; i < number_of_quadrature_points ; i++)
	// 	{	

	// 		// Generate the coordinates
	// 		for  (  int k = 0 ; k < dim ; k++)
	// 		{
	// 			point_coords[k] = quad_triangles->QUAD_POINTS->me[i][k];
	// 		}

	// 		// Create the material points
	// 		//MPS[i] = create_material_point(point_coords,quad_triangles->VOLUMES->ve[i]);

	// 		MPS[i]->beta = beta;
	// 		MPS[i]->kernel_support=RADIAL;
	// 		MPS[i]->neighbours = iv_get(50);


	// 		// Support parameters 
	// 		MPS[i]->MI = m_get(dim,dim);
	// 		MPS[i]->invMI = m_get(dim,dim);

	// 		// set the initial domain of the material points 
	// 		setDomainMaterialPoint(mfree->nodes, grid, MPS[i]);



	// 		MPS[i]->shape_function = NULL;
	// 		MPS[i]->shape_function = mls_shapefunction_materialpoint(MPS[i],2,mfree->nodes);


	// 		// Strain Displacement relationship 
	// 		MPS[i]->B = BMAT(MNULL,MPS[i]->shape_function,dim,IS_AXI,MPS[i]->coords_n_1[0]);


	// 		// Matricies used in internal force calculation
	// 		MPS[i]->fInt = v_get(dim*MPS[i]->num_neighbours);

	// 		// Deformation gradient 
	// 		MPS[i]->F_n = m_get(dim_s,dim_s);
	// 		MPS[i]->inc_F = m_get(dim_s,dim_s);
	// 		m_ident(MPS[i]->inc_F);
	// 		m_ident(MPS[i]->F_n);
	// 		MPS[i]->Jn = 1.00;

	// 		// Temperature
	// 		MPS[i]->temperature = 0;
	// 		if ( mfree->temperatures != NULL){
	// 			for ( int k = 0 ; k < MPS[i]->num_neighbours ; k++)
	// 			{
	// 				int index = MPS[i]->neighbours->ive[k];
	// 				MPS[i]->temperature += MPS[i]->shape_function->phi->ve[k]*mfree->temperatures[index];
	// 			}
	// 		}


	// 		// Density
	// 		MPS[i]->rho = rho;

	// 		// Stress Voigt
	// 		MPS[i]->stressVoigt = v_get(dim_p);

	// 		MPS[i]->temp = m_get(dim_s,dim_s);
	// 		MPS[i]->temp_1 = m_get(dim_s,dim_s);

	// 		// Modify integration factor if problem is axisymmetric
	// 		if ( IS_AXI == 1)
	// 		{
	// 			MPS[i]->INTEGRATION_FACTOR = 2*PI*MPS[i]->coords_n_1[0];
	// 		}else{
	// 			MPS[i]->INTEGRATION_FACTOR = 1.00;
	// 		}


	// 		// Create material state for eahc material point 
	// 		MPS[i]->stateOld= new_material_state(MPS[i]->temperature, material_type,
	// 			dim, IS_AXI);
	// 		MPS[i]->stateNew = new_material_state(MPS[i]->temperature, material_type,
	// 			dim, IS_AXI);


	// 	}



	// 	material_points->MP = MPS;
	// 	material_points->dim = dim;
	// 	material_points->IS_AXI = IS_AXI;





	// }











}
