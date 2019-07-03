#ifndef INTERNAL_FORCE_NEO_H_
#define INTERNAL_FORCE_NEO_H_



#include "internal_force_mooney.h"
#include "Deformation/velocity_grad.h"
#include "Force/Internal/internalForce_Inelastic_Buckley.h"
#include "Deformation/update_Polar_Decomposition.h"
#include "Deformation/rotate_tensor.h"
#include "Material/Buckley/buckleyStress.h"



void internal_force_neo(void *threadarg);

#endif