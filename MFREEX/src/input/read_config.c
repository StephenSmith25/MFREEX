#include "input/read_config.h"












static char *trim(char *str)
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







// Displacement bcs are read into nodesets
static int read_displacementBCs(DOMAIN * domain, FILE * fp)
{
	// INITIALISE 
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * token;
    read = getline(&line, &len, fp);
    line = trim(line);


    // FILE FORMAT 

    // $DisplacementBC
    // name
    // ID
    // disp flags - ( where bcs are applied)
    // disp prescribed ( what value are prescribed)
    // $EndDisplacementBC

    // $
    char * name;
    int ID;
    int BC;

    NODESET * nodeset;
    SIDESET * sideset;
    BLOCKSET * blockset;

	while (strcmp (line, "$ENDDISPLACEMENTBCS") != 0){
      // split string into tokens seperated by equals sign 
      
      while ((strcmp (line, "$EndDisplacementBC") != 0) ){

		  // split string into tokens seperated by equals sign 
	      char * p = strtok(line,"=");

	      // get nodeset by ID
	      if ( strcmp(p, "ID") == 0)
	      {
	      	p = strtok(NULL, " ");
	      	ID = atoi(p);
	      	nodeset = FindNodeSetByID(domain, ID);
	      }


	      // CONSTRAINT IN X DIRECTION
	      if ( strcmp(p, "disp_flag1") == 0)
	      {
	      	p = strtok(NULL, " ");
	      	BC= atoi(p);
	      	if ( BC == 1)
	      	{
	      		//  add a x constraint to the problem
	      		DOF_CONSTRAINT  * dof_constraint = malloc(1*sizeof(DOF_CONSTRAINT));
	      		dof_constraint->next = NULL;

	      		dof_constraint->dir=X;
	      		// add dof constraint to nodeset
	      		AddDOFConstraintToNodeSet(nodeset, dof_constraint);


	      	}
	      }

	      if (( strcmp(p, "disp_ux") == 0 )  && (FindDOFConstraintInNodeSet(nodeset,X) != NULL))
	      {
	      	p = strtok(NULL, " ");
	      	int BC_TYPE= atoi(p);
	      	DOF_CONSTRAINT * dof_constraint = FindDOFConstraintInNodeSet(nodeset,X);
	      	// add dof constraint to nodeset
		    if ( BC_TYPE == 0 )
		    {
		      	dof_constraint->type = DOF_FIXED;
		    }

	      }


	    // CONSTRAINT IN Y DIRECTION
      	if ( strcmp(p, "disp_flag2") == 0)
	      {
	      	p = strtok(NULL, " ");
	      	BC= atoi(p);
	      	if ( BC == 1)
	      	{
	      		//  add a x constraint to the problem
	      		DOF_CONSTRAINT  * dof_constraint = malloc(1*sizeof(DOF_CONSTRAINT));
	      		dof_constraint->next = NULL;
	      		dof_constraint->dir=Y;
	      		// add dof constraint to nodeset
	      		AddDOFConstraintToNodeSet(nodeset, dof_constraint);

	      	}
	      }
	      // check that a y constraint exists
	      if (( strcmp(p, "disp_uy") == 0 )  && (FindDOFConstraintInNodeSet(nodeset,Y) != NULL))
	      {
	      	p = strtok(NULL, " ");
	      	int BC_TYPE= atoi(p);
	      	DOF_CONSTRAINT * dof_constraint = FindDOFConstraintInNodeSet(nodeset,Y);
	      	// add dof constraint to nodeset
		    if ( BC_TYPE == 0 )
		     {
		      	dof_constraint->type = DOF_FIXED;
		    }

	      }

	      // CONSTRAINT IN Z DIRECTION
	      if ( strcmp(p, "disp_flag3") == 0)
	      	{
	      	p = strtok(NULL, " ");
	      	BC= atoi(p);
	      	if ( BC == 1)
	      	{
	      		//  add a x constraint to the problem
	      		DOF_CONSTRAINT  * dof_constraint = malloc(1*sizeof(DOF_CONSTRAINT));
	      		dof_constraint->next = NULL;
	      		dof_constraint->dir=Z;
	      		// add dof constraint to nodeset
	      		AddDOFConstraintToNodeSet(nodeset, dof_constraint);

	      	}
	      }
	      // check that a z constraint exists
	      if (( strcmp(p, "disp_uz") == 0 ) && (FindDOFConstraintInNodeSet(nodeset,Z) != NULL))
	      {
	      	p = strtok(NULL, " ");
	      	int BC_TYPE= atoi(p);
	      	DOF_CONSTRAINT * dof_constraint = FindDOFConstraintInNodeSet(nodeset,Z);
	      	// add dof constraint to nodeset
		    if ( BC_TYPE == 0 )
		     {
		      	dof_constraint->type = DOF_FIXED;
		    }

	      }
		  /* Move to next line  */
	      read = getline(&line, &len, fp);
	      line = trim(line);
 
	}
	  // read next line 
      read = getline(&line, &len, fp);
      line = trim(line);
 
    }

    // FREE ANY ALLOCATED MEMORY
    if (line)
        free(line);

    return 0;



}






int read_config_file(DOMAIN * domain, char * filename)
{


	// initialise 
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * token;
    // mesh_file_name
    char * mesh_file_name;
	// Loop over each  line in the config file 
	// read nodes
    fp = fopen(filename, "r");
    if (fp == NULL){
        fprintf(stderr, "Cannot find config file\n");
        exit(EXIT_FAILURE);
    }
    // LOOP OVER INPUT FILE 
	while ((read = getline(&line, &len, fp)) != -1) {
        line = trim(line);

        if ( strcmp (line, "$DISPLACEMENTBCS")== 0){
        	read_displacementBCs(domain, fp);
        }

        if ( strcmp (line, "$MATERIALS")== 0){
        	//read_nodes(fp, domain);

        	// READ EACH MATERIAL 

        	// READ ID

        	// ADD TO BLOCKSET
        }

        if ( strcmp (line, "$Elements")== 0){
        	//read_elements(fp, domain);
        }
    }




    if (line)
        free(line);

	return 0;

}