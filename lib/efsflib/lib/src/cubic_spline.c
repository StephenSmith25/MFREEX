#include "cubic_spline.h"


int cubic_spline(VEC * w_arr, double r)
{


    if (r <= 0.50 ){
        
        w_arr->ve[0] = 2.0/3 - 4*r*r + 4*r*r*r;
        w_arr->ve[1] = -8*r + 12*r*r;
    }
    else if ((r > 0.50) & (r <= 1.00)){
        w_arr->ve[0] = 4.0/3 - 4*r + 4*r*r - (4.0/3)*r*r*r;
        w_arr->ve[1] = -4 + 8*r - 4*r*r;
    }
    else{        
        w_arr->ve[0] = 0 ;
        w_arr->ve[1] = 0 ;
        
    }
        
    return 0 ;
}

