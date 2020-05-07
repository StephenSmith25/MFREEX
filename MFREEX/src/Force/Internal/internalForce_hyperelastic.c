#include "Force/Internal/internalForce_hyperelastic.h"
static int call_count;
static int print_count;

double internalForce_hyperelastic(VEC *Fint, MATERIAL_POINT *MP, VEC *disp,
                                  VEC *velocity, VEC *matParams,
                                  int (*mat_func_ptr)(VEC *, MAT *, VEC *),
                                  double t_n_1) {

  // Material point
  state_variables *stateNew = MP->stateNew;

  // Get the strain and material parameters from the material point
  double delta_t_min = 1000;
  MAT *F = stateNew->F;
  MAT *Fdot = stateNew->Fdot;
  VEC *stressVoigt = MP->stressVoigt;
  MAT *B;
  IVEC *neighbours;
  MAT *F_r;
  VEC *fInt;

  int num_neighbours = MP->num_neighbours;

  /*  Find deformation gradient */
  B = MP->B;
  neighbours = MP->neighbours;
  F_r = MP->F_n;
  fInt = MP->fInt;
  defgrad_m(MP->inc_F, MP->B, MP->neighbours, num_neighbours, disp);

  // Compute total deformation gradient
  m_mlt(MP->inc_F, MP->F_n, F);

  get_dot_defgrad(Fdot, B, neighbours, F_r, velocity);

  /* Integration parameter */
  double intFactor = MP->volume;
  // Get stress
  mat_func_ptr(stressVoigt, F, matParams);

  // Assemble the force vector
  __zero__(fInt->ve, fInt->max_dim);
  vm_mlt(B, stressVoigt, fInt);

  for (int k = 0; k < num_neighbours; k++) {
    Fint->ve[2 * neighbours->ive[k]] += intFactor * fInt->ve[2 * k];
    Fint->ve[2 * neighbours->ive[k] + 1] += intFactor * fInt->ve[2 * k + 1];
  }

  return delta_t_min;
}
