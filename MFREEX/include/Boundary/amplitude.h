#ifndef AMPLITUDE_H
#define AMPLITUDE_H


#include <stdlib.h>

#include <stdio.h>

typedef enum AMPLITUDE
{
	SMOOTH=1,
	TABULAR=2,
	INSTANT=3
}AMPLITUDE;



typedef struct TABLE
{
	double ** table;

	// x y values in columns 

}TABLE;

typedef struct SMOOTH_AMPLTIUDE
{
	double x_max;
	double x_min;
	double y_max;
	double y_min;

}SMOOTH_AMPLTIUDE;


#endif