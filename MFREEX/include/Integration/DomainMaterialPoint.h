/**********************************
 * @author      Stephen Smith 
 * License:     Not yet decided 
 *
 **********************************/


#ifndef DOMAINMATERIALPOINT_H_
#define DOMAINMATERIALPOINT_H_

#include "mls_shapefunction.h"
#include "Integration/material_point.h"





int setDomainMaterialPoint(meshfreeDomain * mfree, MATERIAL_POINT * MP);

int updateMaterialPointDomain(MATERIAL_POINT * MP);


int neighboursMaterialPoint(MATERIAL_POINT * MP, MAT * nodes);



#endif 