#include "input/read_config.h"



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

        if ( strcmp (line, "$PhysicalNames")== 0){
        	//read_physical_names(fp, domain);
        }

        if ( strcmp (line, "$Nodes")== 0){
        	//read_nodes(fp, domain);
        }

        if ( strcmp (line, "$Elements")== 0){
        	//read_elements(fp, domain);
        }
    }




    if (line)
        free(line);

	return 0;

}