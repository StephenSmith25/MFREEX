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
    
    // READ NEXT LINE
    // FIRST LINE WILL BE FIRST PHYSICAL GROUP WITH FORM
    // {DIM, ID, NAME}
    read = getline(&line, &len, fp);
    line = trim(line);
    int ID = -1;


    // CREATE EACH NODESET
 	domain->NUM_BLOCK_SETS = 0;
 	domain->NUM_SIDE_SETS = 0;
 	domain->NUM_NODE_SETS = 0;
 	domain->sidesets = NULL;
 	domain->nodesets = NULL;
 	domain->blocksets = NULL;


    // LOOP OVER UNTIL END OF PHYSCIAL NAMES 
	while (strcmp (line, "$EndPhysicalNames") != 0){
      char * p = strtok(line," ");
      count = 0;
      while ( p != NULL) 
      {
      	if ( count == 1 )
      		ID = atoi(p);
      	if ( count == 2 ){
      		if ( strncmp(p,"\"PRESSURE\"",4) == 0)
      		{ 
      			// Note: This may be redundent
      			AddPhysicalNameToDomain(domain, PRESSURE_TYPE, ID);

      			// add SIDESET to domain list
      			AddSideSetToDomain(domain, ID);

      		}else if ( strncmp(p,"\"MATERIAL\"",4) == 0)
			{
				// Note: This may be redundent
      			AddPhysicalNameToDomain(domain, PHYSICAL_MATERIAL_TYPE, ID);

      			// add BLOCKSET to domain list
      			AddBlockSetToDomain(domain,ID);
			}else if ( strncmp(p,"\"DISPLACEMENT\"",4) == 0 ){
				//  Note: This may be redundent
      			AddPhysicalNameToDomain(domain, DISPLACEMENT_TYPE, ID);
      			// add NODESET to domain list 
      			AddNodeSetToDomain(domain, ID);			
      		}else{

      		}

      	}
		p = strtok(NULL, " ");

		// INCREMENT COUNTER 
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
    int eDim = 0;
    int ID = -1;

    // pointer to physcial name struct 
    PHYSICAL_NAME * physical_name = NULL;
    BLOCKSET * blockset = NULL;
    NODESET * nodeset = NULL;
    SIDESET * sideset = NULL;
	

	// verticies
	int elementVerticies[20];





    // Loop over each line of elements
	while (strcmp (line, "$EndElements") != 0){

      char * p = strtok(line," ");
      count = 0;
      while ( p != NULL) 
      {

        /* Get element type */ 
      	if ( count == 1 )
      	{
      		eType = atoi(p);
      	}

      	// CHECK WETHER IT IS A NODESET, SIDESET OF MATERIALSET 
      	// fourth entry gives physical type element belongs to 
      	if ( count == 3 ){
      		ID = atoi(p);

      	// FIND physcial_name matching this ID
  			physical_name = FindPhysicalNameByID(domain,ID);

  			// This error should never be reached 
  			if ( physical_name == NULL)
  			{
  				fprintf(stderr,"Cannot find the physical name with ID specified \n");
  				return -1;
  			}

      	}

      	// >6th entry gives element nodes
      	if ( count >= 5 )
      	{


          // if blockset
	      	if ( physical_name->type == PHYSICAL_MATERIAL_TYPE)
	      	{
	      		// Find Blockset that matches this ID;
	      		blockset = FindBlockSetByID(domain, ID);
            
            // create element
             count = 0;
             while ( p != NULL) 
             {
              elementVerticies[count] = atoi(p)-1;

              count++;
              p = strtok(NULL, " ");
             }
             blockset->elements = AddElementToBlockSet(blockset, eType, elementVerticies);


	      	}else if ( physical_name->type == PRESSURE_TYPE)
          {
            sideset = FindSideSetByID(domain, ID);

            // create element
             count = 0;
             while ( p != NULL) 
             {
              elementVerticies[count] = atoi(p)-1;

              count++;
              p = strtok(NULL, " ");
             }
             sideset->elements = AddElementToSideSet(sideset, eType, elementVerticies);

          }else if ( physical_name->type == DISPLACEMENT_TYPE)
          {
            nodeset = FindNodeSetByID(domain, ID);
            // create element
             count = 0;
             while ( p != NULL) 
             {
              elementVerticies[count] = atoi(p)-1;
              count++;
              p = strtok(NULL, " ");
             }
             nodeset->elements = AddElementToNodeSet(nodeset, eType, elementVerticies);

          }else{
              fprintf(stderr,"physical type not supported\n");
          }

          // if sideset 






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



	return domain;
}

