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

	/* TO DO 

	- HANDLE PRECRIBED DISPLACEMENTS
	- HANDLE VELOCITY CONDITIONS.

	*/





	// INITIALISE 
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * token;
    read = getline(&line, &len, fp);
    line = trim(line);



    char * name;
    int ID;
    int BC;

    NODESET * nodeset;

	while (strcmp (line, "$ENDDISPLACEMENTBCS") != 0){
      
      if ((strcmp (line, "$DisplacementBC") == 0) )
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

// material parameters are read into blockset 
static int read_material_parameters(DOMAIN * domain, FILE * fp)
{
	/* TO DO 

	- HANDLE DIRECT APPLICATION
	- HANDLE TABULAR FORM

	*/

	// INITIALISE 
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * token;
    read = getline(&line, &len, fp);
    line = trim(line);
    int num_parameters = 0;
    int ID;

    BLOCKSET * blockset; 

	while (strcmp (line, "$ENDMATERIALS") != 0)
	{
		MATERIAL * new_material = malloc(1*sizeof(MATERIAL));
		if ((strcmp (line, "$Material") == 0) )
		{
      		read = getline(&line, &len, fp);
	      	line = trim(line);

		      while ((strcmp (line, "$EndMaterial") != 0) ){
				
				char * p = strtok(line,"=");

		      	// get blockset by ID
		      	if ( strcmp(p, "ID") == 0)
		      	{
		      		p = strtok(NULL, " ");
		      		ID = atoi(p);
		      		blockset = FindBlockSetByID(domain, ID);
		      	}

		      	if (strcmp(p,"material_name") == 0 )
		      	{


		      		p = strtok(NULL, " ");


		      		if ( strcmp(p,"MAT_RIVLIN") == 0)
		      		{
		      			new_material->material = MAT_RIVLIN;



		      			//new_material->params = v_get(NUM_PARAMETERS);
		      			// next lines should contain the constants
		      			// c1 c2 and kappa

						read = getline(&line, &len, fp);
	      	 	   		line = trim(line);
	      	 	   		 p = strtok(line,"=");

	      	 	   		 // FInd num of paramters

	      	 	   		 if ( strcmp(p,"NUM_PARAMETERS") == 0)
	      	 	   		 {
		      				p = strtok(NULL, " ");
		      				num_parameters = atoi(p);
		      				new_material->params = v_get(num_parameters);

	      	 	   		 }


						read = getline(&line, &len, fp);
	      	 	   		line = trim(line);

	      	 	   		if ( strcmp(line,"$Constants") == 0){
							read = getline(&line, &len, fp);
		      	 	   		line = trim(line);


		      	 	   		while ( strcmp(line,"$EndConstants") != 0)
		      	 	   		{

			      	 	   		p = strtok(line,"=");

			      	 	   		if ( strcmp(p, "C1") == 0)
			      	 	   		{
				      				p = strtok(NULL, " ");
				      				new_material->params->ve[0] = atof(p);
				      			
			      	 	   		}
		 	 	   				 if ( strcmp(p, "C2") == 0)
			      	 	   		 {
									p = strtok(NULL, "");
				      				new_material->params->ve[1] = atof(p);
			      	 	   		 }
			      	 	   		 if ( strcmp(p, "kappa") == 0)
			      	 	   		 {
									p = strtok(NULL, "");
				      				new_material->params->ve[2] = atof(p);
			      	 	   		 }
			      	 	   		read = getline(&line, &len, fp);
			      	 	   		line = trim(line);
		      	 	   		}
	      	 	   		}

	      	 	   		 if ( strcmp(p,"density") == 0)
	      	 	   		 {
		      				p = strtok(NULL, " ");
		      				new_material->density = atof(p);

	      	 	   		 }



						read = getline(&line, &len, fp);
	      	 	   		line = trim(line);


		      		}

		      	}



	      		read = getline(&line, &len, fp);
	      		line = trim(line);

		      }

		}

		blockset->material = new_material;

     	read = getline(&line, &len, fp);
      	line = trim(line);


	}


    // FREE ANY ALLOCATED MEMORY
    if (line)
        free(line);


	return 0;
}

// timestep parameters are read into the timestep structure
static int read_timestep_parameters()
{



	return 0;
}





// pressure conditions are read into sidesets
static int read_pressure_loads(DOMAIN * domain, FILE * fp)
{


	/* TO DO 

	- HANDLE DIRECT APPLICATION
	- HANDLE TABULAR FORM

	*/


	// INITIALISE 
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * token;
    read = getline(&line, &len, fp);
    line = trim(line);

    int ID;


    SIDESET * sideset;

	while (strcmp (line, "$ENDPRESSURES") != 0)
	{
		if ( strcmp(line, "$Pressure") == 0 )
		{
			PRESSURE_LOAD * pressure_load  = malloc(1*sizeof(PRESSURE_LOAD));

	     	while ((strcmp (line, "$EndPressure") != 0) ){



	 		  // split string into tokens seperated by equals sign 
		      char * p = strtok(line,"=");
		     
		      // get nodeset by ID
		      if ( strcmp(p, "ID") == 0)
		      {
		      	p = strtok(NULL, "");
		      	ID = atoi(p);
		      	sideset = FindSideSetByID(domain, ID);
		      	sideset->pressure_load = pressure_load;

		      }
		      if ( strcmp(p, "pressure_type") == 0)
		      {
		      	p = strtok(NULL, "");

		      	if ( strcmp(p,"DIRECT") == 0)
		      	{
		      		sideset->pressure_load->type = DIRECT_PRESSURE;

		      		// move to next line, which contaisn the magnitude
	      	 		read = getline(&line, &len, fp);
	      	 	    line = trim(line);

	      	 	    p = strtok(line,"=");

	      	 	    if ( strcmp(p,"pressure_magnitude") == 0 )
	      	 	    {
	      	 		    p = strtok(NULL, "");
		      			sideset->pressure_load->pMag = atof(p);
	      	 	    }

	      	 	    // move to next line, which contains the ramp
	      	 		read = getline(&line, &len, fp);
	      	 	    line = trim(line);

	      	 	    p = strtok(line,"=");
  	 	    		if ( strcmp(p,"pressure_ramp") == 0)
	      	 	    {

	      	 		    p = strtok(NULL, "");
	      	 		    if ( strcmp(p,"SMOOTH") == 0)
	      	 		    {
	      	 		    	// Set up smooth table 
	      	 		    	sideset->pressure_load->amplitude=SMOOTH;
	      	 		    	sideset->pressure_load->smooth_table = malloc(1*sizeof(SMOOTH_AMPLTIUDE));
	      	 		    	
	      	 		    	/* Read smooth points i.e x0,y0 x1,y1 into load struct */
  	 						read = getline(&line, &len, fp);
	      	 	   			read = getline(&line, &len, fp);
	      	 	   			line = trim(line);
	      	 	   			p = strtok(line," ");
	      	 	   			sideset->pressure_load->smooth_table->x_min = atof(p);
	      	 		  		p = strtok(NULL, " ");
	      	 	   			sideset->pressure_load->smooth_table->y_min = atof(p);
							read = getline(&line, &len, fp);
	      	 	   			line = trim(line);
	      	 	   			p = strtok(line," ");
	      	 	   			sideset->pressure_load->smooth_table->x_max = atof(p);
	      	 		  		p = strtok(NULL, " ");
	      	 	   			sideset->pressure_load->smooth_table->y_max = atof(p);
	      	 		    }


	      	 	    }




		      	}

		      }

		      // read next line 
	      	  read = getline(&line, &len, fp);
	      	  line = trim(line);

			}
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
	        read_material_parameters(domain, fp);
        }

        if ( strcmp (line, "$PRESSURES")== 0){
        	read_pressure_loads(domain,fp);
      	 }
    }




    if (line)
        free(line);

	return 0;

}