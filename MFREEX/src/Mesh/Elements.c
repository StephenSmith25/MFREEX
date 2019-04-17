#include "Mesh/Elements.h"


// Create a new element 
ELEMENT * CreateNewElement(ELEMENT_TYPE etype, int * verticies)
{
	ELEMENT * element = malloc(1*sizeof(ELEMENT));
	element->etype = etype;


	if ( etype == TRIANGLES)
	{
		element->verticies = malloc(3*sizeof(int));
		element->verticies[0] = verticies[0];
		element->verticies[1] = verticies[1];
		element->verticies[2] = verticies[2];

	}else if ( etype == LINE)
	{
		element->verticies = malloc(3*sizeof(int));

		element->verticies[0] = verticies[0];
		element->verticies[1] = verticies[1];

	}else if ( etype == QUAD)
	{
		element->verticies = malloc(4*sizeof(int));

		element->verticies[0] = verticies[0];
		element->verticies[1] = verticies[1];
		element->verticies[2] = verticies[2];
		element->verticies[3] = verticies[3];

	}else if ( etype == TETRAHEDRON)
	{
		element->verticies = malloc(4*sizeof(int));

		element->verticies[0] = verticies[0];
		element->verticies[1] = verticies[1];
		element->verticies[2] = verticies[2];
		element->verticies[3] = verticies[3];


	}else{
		fprintf(stderr, "ELEMENT not supported\n");
	}


	return element;


}

void DeleteElement(ELEMENT * element)
{
	free(element);
}