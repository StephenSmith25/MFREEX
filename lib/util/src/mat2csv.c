#include "mat2csv.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>



void mat2csv ( MAT * in ,char * folderName, const char * fileName)

{

	FILE * fp ;


	struct stat st = {0};

	if ( stat(folderName, &st) == -1) {
		mkdir(folderName, 07777);
	}



	char buffer[50] ;
	snprintf(buffer, 50, "%s%s%s",folderName,"/",fileName);
	
	fp = fopen(buffer,"w");


	for( int i = 0 ; i < in->m ;  i++){

		for(int j = 0 ; j < in->n ;  j++){
			fprintf(fp,"%lf",in->me[i][j]);


			if ( j < in->n-1){
				fprintf(fp,",");
			}



		}


		fprintf(fp,"\n");


	}
	



	fclose(fp);
	return ; 









}

void disp2csv ( MAT * in ,char * folderName, const char * fileName)

{

	FILE * fp ;


	struct stat st = {0};

	if ( stat(folderName, &st) == -1) {
		mkdir(folderName, 07777);
	}
	


	char buffer[50] ;
	snprintf(buffer, 50, "%s%s%s",folderName,"/",fileName);
	
	fp = fopen(buffer,"w");


	for( int i = 0 ; i < in->m ;  i++){

		if ( i == 0)
		fprintf(fp,"x,y\n");

		for(int j = 0 ; j < in->n ;  j++){
			fprintf(fp,"%lf",in->me[i][j]);


			if ( j < in->n-1){
				fprintf(fp,",");
			}



		}


		fprintf(fp,"\n");


	}
	



	fclose(fp);
	return ; 









}
