
 #include <sys/stat.h>
 #include <unistd.h>


#include "PostProcess/saveDisp.h"

void saveDisp(MAT * in, state_variables ** state, char * folderName, const char * fileName)

{

	FILE * fp ;


	struct stat st = {0};

	if ( stat(folderName, &st) == -1) {
		mkdir(folderName, 07777);
	}
	


	char buffer[50] ;
	snprintf(buffer, 50, "%s%s%s",folderName,"/",fileName);
	
	fp = fopen(buffer,"w");
	double sx,sy,sz =0;


	for( int i = 0 ; i < in->m ;  i++){


		MAT * sigma = state[i]->sigma;

		if ( i == 0){
		fprintf(fp,"r,z,sr,sz,s0\n");
		}


		for(int j = 0 ; j < in->n ;  j++){
			fprintf(fp,"%lf",in->me[i][j]);
			if ( j < in->n){
				fprintf(fp,",");
			}
		}

		for(int j = 0 ; j < 3 ;  j++){
			fprintf(fp,"%lf",sigma->me[j][j]);


			if ( j < 3){
				fprintf(fp,",");
			}
		}


		fprintf(fp,"\n");


	}
	



	fclose(fp);
	return ; 









}