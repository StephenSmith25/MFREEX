#include "Mesh/Elements.h"


// Create a new element 
ELEMENT * CreateNewElement(ELEMENT_TYPE etype, double * verticies)
{
	ELEMENT * element = malloc(1*sizeof(ELEMENT));
	element->etype = etype;
	return element;


}