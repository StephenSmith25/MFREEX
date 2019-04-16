#ifndef ELEMENTS_H_
#define ELEMENTS_H_



#include <stdlib.h>
#include <stdio.h>



// Element types
typedef enum ELEMENT_TYPE
{
	LINE = 1,
	TRIANGLES=2,
	QUAD=3,
	TETRAHEDRON=4
}ELEMENT_TYPE; 

// ELEMENT OBJECT
typedef struct ELEMENT
{
	ELEMENT_TYPE etype;
	int * verticies;

}ELEMENT;

typedef struct ELEMENT_LIST
{
	ELEMENT * element;

}ELEMENT_LIST;


#endif 