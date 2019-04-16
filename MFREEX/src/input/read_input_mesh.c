#include "input/read_input_mesh.h"


// from https://stackoverflow.com/questions/122616/how-do-i-trim-leading-trailing-whitespace-in-a-standard-way?page=1&tab=votes#tab-top
// credit https://stackoverflow.com/users/19719/indiv
char *trim(char *str)
{
    size_t len = 0;
    char *frontp = str;
    char *endp = NULL;

    if( str == NULL ) { return NULL; }
    if( str[0] == '\0' ) { return str; }

    len = strlen(str);
    endp = str + len;

    /* Move the front and back pointers to address the first non-whitespace
     * characters from each end.
     */
    while( isspace((unsigned char) *frontp) ) { ++frontp; }
    if( endp != frontp )
    {
        while( isspace((unsigned char) *(--endp)) && endp != frontp ) {}
    }

if( str + len - 1 != endp )
    *(endp + 1) = '\0';
else if( frontp != str &&  endp == frontp )
    *str = '\0';

    /* Shift the string so that it starts at str so that if it's dynamically
     * allocated, we can still free it on the returned pointer.  Note the reuse
     * of endp to mean the front of the string buffer now.
     */
endp = str;
if( frontp != str )
{
    while( *frontp ) { *endp++ = *frontp++; }
    *endp = '\0';
}


return str;
}

static inline void read_physical_names(FILE * fp, DOMAIN * domain )
{

	// INITIALISE AND GET NUMBER OF PHYSICAL GROUPS 
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * token;
    int line_count = 0;
    int count = 0;
    read = getline(&line, &len, fp);
    line = trim(line);
    int number_of_physical_groups = atoi(line);
    domain->number_of_physical_groups = number_of_physical_groups;
    PHYSICAL_GROUP * physical_groups = malloc(number_of_physical_groups * sizeof(PHYSICAL_GROUP));
    

    printf("number of physcial groups = %d \n", number_of_physical_groups);
    // READ NEXT LINE
    read = getline(&line, &len, fp);
    line = trim(line);
    int ID = -1;
    // LOOP OVER UNTIL END OF PHYSCIAL NAMES 
	while (strcmp (line, "$EndPhysicalNames") != 0){
      char * p = strtok(line," ");
      count = 0;
      while ( p != NULL) 
      {
      	if ( count == 1 )
      		ID = atoi(p);

      	if ( count == 2 ){
      		if ( strcmp(p,"PRESSURE") == 0)
      		{
      			physical_groups[ID-1] = PRESSURE_TYPE;
      		}else if (strcmp(p,"PRESSURE") == 0)
			{
      			physical_groups[ID-1] = PHYSICAL_MATERIAL_TYPE;
			}else{
      			physical_groups[ID-1] = DISPLACEMENT_TYPE;
			}

      	}
		p = strtok(NULL, " ");


		count++;
      }

      read = getline(&line, &len, fp);
      line = trim(line);
      ++line_count;
 
    }

    // FREE ANY ALLOCATED MEMORY
    if (line)
        free(line);

    return ;

}


int read_elements(FILE * fp, DOMAIN * domain)
{
	// INITIALISE AND GET NUMBER OF PHYSICAL GROUPS 
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * token;
    int line_count = 0;
    int count = 0;
    read = getline(&line, &len, fp);
    line = trim(line);

    int total_number_of_elements = atoi(line);
    read = getline(&line, &len, fp);
    line = trim(line);

    int eType = 0;
    int pType = 0;

    // Loop over each line of elements
	while (strcmp (line, "$EndElements") != 0){

      char * p = strtok(line," ");
      count = 0;
      while ( p != NULL) 
      {

      	// second entry gives element type
      	// 1 = line
      	// 2 = triangle
      	// 3 = quad
      	// 4 = tet 
      	if ( count == 1 )
      	{
      		eType = atoi(p);
      	}

      	// fourth entry gives physical type element belongs to 
      	if ( count == 3 ){
      		pType = atoi(p);
      	}

      	if ( pType == PHYSICAL_MATERIAL_TYPE)
      	{

      	}
      	// >6th entry gives element nodes
      	if ( count > 5 )
      	{

      	}


      	// move to next token in string
		p = strtok(NULL, " ");
		count++;
      }

      read = getline(&line, &len, fp);
      line = trim(line);
      ++line_count;
 
    }

	// loop over each physcial group and find the elements belonging to each group
	
	int ID = 0;
	// FIND WHICH GROUP THESE NODES BELONG TO 
	// 5TH ENTRY IN EACH LINE CONTAINS THE ID WHICH REPRESENTS THE PHYSCIAL GROUP OF (ID)

	if ( ID == PHYSICAL_MATERIAL_TYPE)
	{
		// create a blockset

		// ELEMENT TYPE 

	}

	// if ( ID == DISPLACEMENT_TYPE)
	// {
	// 	// create a nodeset

	// 	// ELEMENT TYPE
	// }

	// if ( ID == PRESSURE_TYPE)
	// {
	// 	// create a pressure set

	// 	// ELEMENT TYPE
	// }

	return 0;
}


int read_nodes(FILE * fp, DOMAIN * domain)
{
	// initialise 
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * token;
    int line_count = 0;
    int count = 0;
    read = getline(&line, &len, fp);
    line = trim(line);
    int number_of_nodes = atoi(line);
    MAT * nodes = m_get(number_of_nodes,DIM);
    
    read = getline(&line, &len, fp);
    line = trim(line);


	while (strcmp (line, "$EndNodes") != 0){
      char * p = strtok(line," ");
      count = 0;
      while ( p != NULL) 
      {
      	if ( count > 0 ){
      		nodes->me[line_count][count-1] = atof(p);
      	}
		p = strtok(NULL, " ");


		count++;
      }

      read = getline(&line, &len, fp);
      line = trim(line);
      ++line_count;
 
    }

    if (line)
        free(line);

   	domain->NODES = nodes;

    return 0;
}


DOMAIN * read_intput_mesh(char * input_file_name)
{


	// initialise 
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * token;
    // mesh_file_name
    char * mesh_file_name;

    MAT * elements;
    MAT * nodes;

	
	DOMAIN * domain = malloc(1*sizeof(DOMAIN));
	// read nodes
    fp = fopen(input_file_name, "r");
    if (fp == NULL){
        fprintf(stderr, "Cannot find input file\n");
        exit(EXIT_FAILURE);
    }
    // LOOP OVER INPUT FILE 
	while ((read = getline(&line, &len, fp)) != -1) {
        line = trim(line);

        if ( strcmp (line, "$PhysicalNames")== 0){
        	read_physical_names(fp, domain);
        }

        if ( strcmp (line, "$Nodes")== 0){
        	read_nodes(fp, domain);
        }

        if ( strcmp (line, "$Elements")== 0){
        	read_elements(fp, domain);
        }
    }




    if (line)
        free(line);


    // read config file to get boundary conditions and 

	return domain;
}

