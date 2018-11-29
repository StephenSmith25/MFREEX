#include "smoothstep.h"
double smoothstep(double x, double x_max, double x_min){


double normal_x = (x - x_min)/(x_max - x_min) ;

if ( normal_x <= 0 ){
	
	return 0 ;
}

else if (( normal_x > 0 ) && (normal_x <= 1)){

	return 3*normal_x*normal_x  - 2* normal_x*normal_x*normal_x ; 
}

else {
return 1; 
}



}
