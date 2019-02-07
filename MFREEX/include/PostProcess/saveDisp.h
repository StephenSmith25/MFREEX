#ifndef SAVEDISP_H_
#define SAVEDISP_H_


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "matrix.h"
#include "matrix2.h"
#include "Material/material.h"






void saveDisp(MAT * in, state_variables ** state,char * folderName, const char * filename);


#endif