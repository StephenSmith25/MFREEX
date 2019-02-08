#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "matrix2.h"
#include <sys/time.h>
#include <omp.h>

const int NUM_TRIES =10;
const int NUM_MATRIX_MULTI = 10;
const int NUM_POINTS = 300;
const int NUM_THREADS = 4;
int do_stuff(MAT * A, MAT*B, MAT * C, VEC * largeVector);

int main(int argc, char const *argv[])
{

	omp_set_num_threads(4);






	struct timeval start, end;
	// generate clipped voronoi diagram
	gettimeofday(&start, NULL);
#pragma omp parallel
{

	int n = 0 ;
	VEC * largeArray = v_get(50);
	MAT * A = m_get(3,3);
	MAT * B = m_get(3,3);
	MAT * C = m_get(3,3);

	while ( n < NUM_TRIES){
		#pragma omp for 
		for ( int i = 0 ; i < NUM_POINTS ; i++)
		{
			if ( i == 1)
			{
				printf("i = 1\n");
			}
			do_stuff(A, B, C, largeArray);
		}
	
	#pragma omp atomic
		++n;


	}




}


 	// get time took to run
	gettimeofday(&end, NULL);
	double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6;

	// get time taken to run
	printf("while loop took %lf seconds to run\n", delta);



	/* code */
	return 0;
}


int do_stuff(MAT * A,MAT * B, MAT * C ,VEC * largeVector )
{

	for ( int i = 0 ; i < NUM_MATRIX_MULTI ; i++)
	{
		m_mlt(A,B,C);
	}



	return 0;
}