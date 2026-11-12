/*
 * Kinematics and constitutive models.
 *
 *  - poldec: polar decomposition F = R U = V R
 *  - SVK:    Saint Venant-Kirchhoff hyperelasticity (returns 1st Piola-
 *            Kirchhoff stress in Voigt order [P11, P22, P12, P21])
 *
 * Checks are physical identities: zero stress at zero strain, analytic
 * uniaxial response, objectivity under rigid rotation.
 */
#include <math.h>

#include "matrix.h"
#include "matrix2.h"

#include "Deformation/poldec.h"

#include "mfreex_test.h"

/* Declared here rather than via its header: the header is renamed on the
 * modernize branch (SVK.h -> svk_defgrad.h), the function is not. */
int SVK(VEC *stressVoigt, MAT *defGrad, VEC *params);

#define TOL 1e-12

/* ------------------------------------------------------------- poldec -- */

static void check_polar(MAT *F) {
  int n = F->m;
  MAT *R = m_get(n, n), *U = m_get(n, n), *V = m_get(n, n);
  MAT *T = m_get(n, n);
  poldec(F, R, U, V);

  m_mlt(R, U, T); /* F = R U */
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) CHECK_CLOSE(T->me[i][j], F->me[i][j], 1e-10);

  mtrm_mlt(R, R, T); /* R^T R = I */
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) CHECK_CLOSE(T->me[i][j], i == j ? 1.0 : 0.0, 1e-10);

  for (int i = 0; i < n; ++i) /* U, V symmetric */
    for (int j = 0; j < n; ++j) {
      CHECK_CLOSE(U->me[i][j], U->me[j][i], 1e-10);
      CHECK_CLOSE(V->me[i][j], V->me[j][i], 1e-10);
    }
}

TEST_CASE(polar_decomposition_2d) {
  MAT *F = m_get(2, 2);
  F->me[0][0] = 1.3; F->me[0][1] = 0.4;
  F->me[1][0] = -0.2; F->me[1][1] = 0.9;
  check_polar(F);
}

TEST_CASE(polar_decomposition_3d) {
  MAT *F = m_get(3, 3);
  const double v[9] = {1.1, 0.2, 0.0, -0.3, 0.95, 0.1, 0.05, 0.0, 1.2};
  for (int i = 0; i < 9; ++i) F->me[i / 3][i % 3] = v[i];
  check_polar(F);
}

TEST_CASE(polar_decomposition_of_pure_rotation) {
  double t = 0.7;
  MAT *F = m_get(2, 2), *R = m_get(2, 2), *U = m_get(2, 2), *V = m_get(2, 2);
  F->me[0][0] = cos(t); F->me[0][1] = -sin(t);
  F->me[1][0] = sin(t); F->me[1][1] = cos(t);
  poldec(F, R, U, V);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      CHECK_CLOSE(U->me[i][j], i == j ? 1.0 : 0.0, 1e-10);
      CHECK_CLOSE(R->me[i][j], F->me[i][j], 1e-10);
    }
}

/* ---------------------------------------------------------------- SVK -- */

static const double LAMBDA = 1200.0, MU = 800.0;

static VEC *svk_params(void) {
  VEC *p = v_get(2);
  p->ve[0] = LAMBDA;
  p->ve[1] = MU;
  return p;
}

static VEC *svk_stress_2d(double F11, double F12, double F21, double F22) {
  MAT *F = m_get(2, 2);
  F->me[0][0] = F11; F->me[0][1] = F12;
  F->me[1][0] = F21; F->me[1][1] = F22;
  VEC *P = v_get(4);
  SVK(P, F, svk_params());
  return P;
}

TEST_CASE(svk_zero_stress_at_identity) {
  VEC *P = svk_stress_2d(1, 0, 0, 1);
  for (int i = 0; i < 4; ++i) CHECK_CLOSE(P->ve[i], 0.0, TOL);
}

TEST_CASE(svk_uniaxial_stretch_matches_analytic) {
  /* F = diag(l, 1): E11 = (l^2-1)/2, S11 = (lambda+2mu) E11,
   * S22 = lambda E11, P = F S */
  double l = 1.25, E11 = 0.5 * (l * l - 1);
  VEC *P = svk_stress_2d(l, 0, 0, 1);
  CHECK_CLOSE(P->ve[0], l * (LAMBDA + 2 * MU) * E11, TOL);
  CHECK_CLOSE(P->ve[1], LAMBDA * E11, TOL);
  CHECK_CLOSE(P->ve[2], 0.0, TOL);
  CHECK_CLOSE(P->ve[3], 0.0, TOL);
}

TEST_CASE(svk_simple_shear_matches_analytic) {
  /* F = [[1, g], [0, 1]]: C = [[1, g], [g, 1+g^2]],
   * E = [[0, g/2], [g/2, g^2/2]] */
  double g = 0.3;
  double E12 = 0.5 * g, E22 = 0.5 * g * g, trE = E22;
  double S11 = LAMBDA * trE, S12 = 2 * MU * E12, S22 = LAMBDA * trE + 2 * MU * E22;
  VEC *P = svk_stress_2d(1, g, 0, 1);
  CHECK_CLOSE(P->ve[0], S11 + g * S12, TOL); /* P11 */
  CHECK_CLOSE(P->ve[1], S22, TOL);           /* P22 */
  CHECK_CLOSE(P->ve[2], S12 + g * S22, TOL); /* P12 */
  CHECK_CLOSE(P->ve[3], S12, TOL);           /* P21 */
}

TEST_CASE(svk_rigid_rotation_is_stress_free) {
  double t = 0.9;
  VEC *P = svk_stress_2d(cos(t), -sin(t), sin(t), cos(t));
  for (int i = 0; i < 4; ++i) CHECK_CLOSE(P->ve[i], 0.0, 1e-9);
}

TEST_CASE(svk_objectivity_P_rotates_with_F) {
  /* For F' = Q F:  P' = Q P */
  double t = 0.4, c = cos(t), s = sin(t);
  double F11 = 1.2, F12 = 0.15, F21 = -0.05, F22 = 0.9;
  VEC *P = svk_stress_2d(F11, F12, F21, F22);
  VEC *Pr = svk_stress_2d(c * F11 - s * F21, c * F12 - s * F22,
                          s * F11 + c * F21, s * F12 + c * F22);
  double P11 = P->ve[0], P22 = P->ve[1], P12 = P->ve[2], P21 = P->ve[3];
  CHECK_CLOSE(Pr->ve[0], c * P11 - s * P21, 1e-9);
  CHECK_CLOSE(Pr->ve[2], c * P12 - s * P22, 1e-9);
  CHECK_CLOSE(Pr->ve[3], s * P11 + c * P21, 1e-9);
  CHECK_CLOSE(Pr->ve[1], s * P12 + c * P22, 1e-9);
}

int main(void) {
  RUN(polar_decomposition_2d);
  RUN(polar_decomposition_3d);
  RUN(polar_decomposition_of_pure_rotation);
  RUN(svk_zero_stress_at_identity);
  RUN(svk_uniaxial_stretch_matches_analytic);
  RUN(svk_simple_shear_matches_analytic);
  RUN(svk_rigid_rotation_is_stress_free);
  RUN(svk_objectivity_P_rotates_with_F);
  return TEST_SUMMARY();
}
