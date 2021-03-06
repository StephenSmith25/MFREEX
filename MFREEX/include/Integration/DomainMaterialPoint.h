/**********************************
 * @author      Stephen Smith 
 * License:     Not yet decided 
 *
 **********************************/


#ifndef DOMAINMATERIALPOINT_H_
#define DOMAINMATERIALPOINT_H_

#include "mls_shapefunction.h"
#include "Integration/material_point.h"
#include "m_inverse_small.h"
#include "cellSearch.h"





int setDomainMaterialPoint(MAT * nodes,  CELLS * cells, MATERIAL_POINT * MP);

int updateDomainMaterialPoint(MAT * nodes,CELLS * cells,  MATERIAL_POINT * MP);


int neighboursMaterialPoint(MATERIAL_POINT * MP, MAT * nodes);



#endif 