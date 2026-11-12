/*
 * Hyperelastic constitutive models -- energy consistency.
 *
 * For a hyperelastic law the first Piola-Kirchhoff stress is P = dW/dF.
 * Each model is checked against a central finite difference of its strain
 * energy W(F), written independently below. This catches wrong factors,
 * transposes and missing terms, and (with P(I) = 0) the reference state.
 *
 * 2D models are plane strain: F is 2x2 with F33 = 1 implied.
 * Axisymmetric/3D use a 3x3 F with no out-of-plane shear; Voigt order is
 * [P11, P22, P12, P21, (P33)].
 */
#include <math.h>

#include "matrix.h"
#include "matrix2.h"

#include "Material/Hyperelastic/compress_neo.h"
#include "Material/Hyperelastic/cubicRivlin.h"
#include "Material/Hyperelastic/mooneyRivlin.h"
#include "Material/Hyperelastic/svk.h"
#include "Material/Hyperelastic/yeoh.h"
#include "Material/material.h"

#include "mfreex_test.h"

/* ---- 3x3 helpers on a full F (2x2 input is embedded with F33 = 1) ---- */
typedef struct { double a[3][3]; } T3;

static T3 embed(MAT *F) {
  T3 t = {{{0}}};
  t.a[2][2] = 1;
  for (int i = 0; i < (int)F->m; ++i)
    for (int j = 0; j < (int)F->n; ++j) t.a[i][j] = F->me[i][j];
  return t;
}
static T3 rcg(T3 F) { /* C = F^T F */
  T3 C = {{{0}}};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) C.a[i][j] += F.a[k][i] * F.a[k][j];
  return C;
}
static double tr(T3 A) { return A.a[0][0] + A.a[1][1] + A.a[2][2]; }
static double tr2(T3 A) { /* tr(A A) */
  double s = 0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) s += A.a[i][j] * A.a[j][i];
  return s;
}
static double det(T3 A) {
  return A.a[0][0] * (A.a[1][1] * A.a[2][2] - A.a[1][2] * A.a[2][1]) -
         A.a[0][1] * (A.a[1][0] * A.a[2][2] - A.a[1][2] * A.a[2][0]) +
         A.a[0][2] * (A.a[1][0] * A.a[2][1] - A.a[1][1] * A.a[2][0]);
}

/* ---- strain energies ---- */
typedef double (*energy_fn)(T3 F, const double *p);

static double W_svk(T3 F, const double *p) { /* p = {lambda, mu} */
  T3 C = rcg(F), E;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) E.a[i][j] = 0.5 * (C.a[i][j] - (i == j));
  return 0.5 * p[0] * tr(E) * tr(E) + p[1] * tr2(E);
}
static double W_neo(T3 F, const double *p) { /* p = {mu, lambda} */
  T3 C = rcg(F);
  double lnJ = log(det(F));
  return 0.5 * p[0] * (tr(C) - 3) - p[0] * lnJ + 0.5 * p[1] * lnJ * lnJ;
}
static double W_mooney(T3 F, const double *p) { /* p = {c1, c2, lambda} */
  T3 C = rcg(F);
  double I1 = tr(C), I2 = 0.5 * (I1 * I1 - tr2(C)), I3 = det(C);
  double p0 = p[0] + 2 * p[1];
  return p[0] * (I1 - 3) + p[1] * (I2 - 3) - p0 * log(I3) + 0.5 * p[2] * log(I3) * log(I3);
}
static double W_iso_poly(T3 F, const double *p) { /* p = {c10, c20, c30, K} */
  T3 C = rcg(F);
  double J = det(F), I1b = tr(C) * pow(J, -2.0 / 3.0), x = I1b - 3;
  return p[0] * x + p[1] * x * x + p[2] * x * x * x + 0.5 * p[3] * (J - 1) * (J - 1);
}

/* ---- calling the models ---- */
typedef enum { STATE_SVK, STATE_NEO, STATE_MOONEY, DEFGRAD_YEOH, DEFGRAD_CUBIC } model;

static void stress(model m, MAT *F, const double *params, int np, VEC *P) {
  VEC *prm = v_get(np);
  for (int i = 0; i < np; ++i) prm->ve[i] = params[i];
  if (m == DEFGRAD_YEOH) { yeoh(P, F, prm); return; }
  if (m == DEFGRAD_CUBIC) { cubicRivlin(P, F, prm); return; }
  int axi = F->m == 3;
  state_variables *s = new_material_state(20.0, HYPERELASTIC, 2, axi);
  m_copy(F, s->F);
  s->Jacobian = det(embed(F));
  if (m == STATE_SVK) svk(P, s, prm);
  if (m == STATE_NEO) compress_neo(P, s, prm);
  if (m == STATE_MOONEY) mooneyRivlin(P, s, prm);
}

static const int VI[5] = {0, 1, 0, 1, 2}, VJ[5] = {0, 1, 1, 0, 2};

static void check_energy(model m, energy_fn W, const double *p, int np, MAT *F) {
  int nv = F->m == 3 ? 5 : 4;
  VEC *P = v_get(nv);
  stress(m, F, p, np, P);
  const double h = 1e-6;
  for (int v = 0; v < nv; ++v) {
    T3 Fp = embed(F), Fm = embed(F);
    Fp.a[VI[v]][VJ[v]] += h;
    Fm.a[VI[v]][VJ[v]] -= h;
    double dW = (W(Fp, p) - W(Fm, p)) / (2 * h);
    CHECK_CLOSE(P->ve[v], dW, 1e-6);
  }
}

static void check_reference_state(model m, const double *p, int np, int n) {
  MAT *F = m_get(n, n);
  m_ident(F);
  VEC *P = v_get(n == 3 ? 5 : 4);
  stress(m, F, p, np, P);
  for (int v = 0; v < (int)P->dim; ++v) CHECK_CLOSE(P->ve[v], 0.0, 1e-10);
}

static MAT *F2(double a, double b, double c, double d) {
  MAT *F = m_get(2, 2);
  F->me[0][0] = a; F->me[0][1] = b; F->me[1][0] = c; F->me[1][1] = d;
  return F;
}
static MAT *F3(double a, double b, double c, double d, double hoop) {
  MAT *F = m_get(3, 3);
  F->me[0][0] = a; F->me[0][1] = b; F->me[1][0] = c; F->me[1][1] = d;
  F->me[2][2] = hoop;
  return F;
}

/* A few deformations: stretch, shear, general (incl. rotation) */
#define FOR_EACH_F2(body) do { MAT *F; \
  F = F2(1.2, 0.0, 0.0, 0.9); body; \
  F = F2(1.0, 0.3, 0.0, 1.0); body; \
  F = F2(1.1, 0.25, -0.15, 0.95); body; } while (0)
#define FOR_EACH_F3(body) do { MAT *F; \
  F = F3(1.2, 0.0, 0.0, 0.9, 1.05); body; \
  F = F3(1.0, 0.3, 0.0, 1.0, 1.0); body; \
  F = F3(1.1, 0.25, -0.15, 0.95, 1.1); body; } while (0)

static const double P_SVK[2] = {1200.0, 800.0};
static const double P_NEO[2] = {800.0, 1200.0};
static const double P_MR[3] = {300.0, 50.0, 1000.0};
static const double P_POLY[4] = {400.0, -20.0, 5.0, 5000.0};

TEST_CASE(svk_plane_strain)        { check_reference_state(STATE_SVK, P_SVK, 2, 2); FOR_EACH_F2(check_energy(STATE_SVK, W_svk, P_SVK, 2, F)); }
TEST_CASE(svk_axisymmetric)        { check_reference_state(STATE_SVK, P_SVK, 2, 3); FOR_EACH_F3(check_energy(STATE_SVK, W_svk, P_SVK, 2, F)); }
TEST_CASE(neo_hookean_plane_strain){ check_reference_state(STATE_NEO, P_NEO, 2, 2); FOR_EACH_F2(check_energy(STATE_NEO, W_neo, P_NEO, 2, F)); }
TEST_CASE(neo_hookean_axisymmetric){ check_reference_state(STATE_NEO, P_NEO, 2, 3); FOR_EACH_F3(check_energy(STATE_NEO, W_neo, P_NEO, 2, F)); }
TEST_CASE(mooney_rivlin_plane_strain){ check_reference_state(STATE_MOONEY, P_MR, 3, 2); FOR_EACH_F2(check_energy(STATE_MOONEY, W_mooney, P_MR, 3, F)); }
TEST_CASE(mooney_rivlin_axisymmetric){ check_reference_state(STATE_MOONEY, P_MR, 3, 3); FOR_EACH_F3(check_energy(STATE_MOONEY, W_mooney, P_MR, 3, F)); }
TEST_CASE(yeoh_3d)                 { check_reference_state(DEFGRAD_YEOH, P_POLY, 4, 3); FOR_EACH_F3(check_energy(DEFGRAD_YEOH, W_iso_poly, P_POLY, 4, F)); }
TEST_CASE(cubic_rivlin_3d)         { check_reference_state(DEFGRAD_CUBIC, P_POLY, 4, 3); FOR_EACH_F3(check_energy(DEFGRAD_CUBIC, W_iso_poly, P_POLY, 4, F)); }
TEST_CASE(cubic_rivlin_plane_strain){ check_reference_state(DEFGRAD_CUBIC, P_POLY, 4, 2); FOR_EACH_F2(check_energy(DEFGRAD_CUBIC, W_iso_poly, P_POLY, 4, F)); }

int main(void) {
  RUN(svk_plane_strain);
  RUN(svk_axisymmetric);
  RUN(neo_hookean_plane_strain);
  RUN(neo_hookean_axisymmetric);
  RUN(mooney_rivlin_plane_strain);
  RUN(mooney_rivlin_axisymmetric);
  RUN(yeoh_3d);
  RUN(cubic_rivlin_3d);
  RUN(cubic_rivlin_plane_strain);
  return TEST_SUMMARY();
}
