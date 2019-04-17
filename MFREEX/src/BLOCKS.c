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
