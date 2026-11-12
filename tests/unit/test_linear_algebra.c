/*
 * Small dense linear algebra used throughout MFREEX:
 *   - util: determinant, trace, m_inverse_small (2x2 / 3x3)
 *   - eig3 (lib/Eigen): Kopp's dsyevh3 / dsyevq3 / dsyevc3 3x3 symmetric
 *     eigen-solvers, used by the Buckley material model and polar
 *     decomposition.
 */
#include "matrix.h"
#include "matrix2.h"

#include "determinant.h"
#include "m_inverse_small.h"
#include "trace.h"
#include "contraction.h"
#include "matop_3x3.h"
#include "smoothstep.h"

#include "dsyevc3.h"
#include "dsyevh3.h"
#include "dsyevq3.h"

#include "mfreex_test.h"

#define TOL 1e-12

static MAT *mat2(double a, double b, double c, double d) {
  MAT *A = m_get(2, 2);
  A->me[0][0] = a; A->me[0][1] = b;
  A->me[1][0] = c; A->me[1][1] = d;
  return A;
}

static MAT *mat3(const double v[9]) {
  MAT *A = m_get(3, 3);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) A->me[i][j] = v[3 * i + j];
  return A;
}

/* ---------------------------------------------------------------- util -- */

TEST_CASE(determinant_2x2_and_3x3) {
  CHECK_CLOSE(determinant(mat2(3, 1, 2, 4)), 10.0, TOL);
  const double a[9] = {2, -3, 1, 2, 0, -1, 1, 4, 5};
  CHECK_CLOSE(determinant(mat3(a)), 49.0, TOL);
  const double sing[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  CHECK_CLOSE(determinant(mat3(sing)), 0.0, TOL);
}

TEST_CASE(trace_3x3) {
  const double a[9] = {2, -3, 1, 2, 7, -1, 1, 4, 5};
  CHECK_CLOSE(trace(mat3(a)), 14.0, TOL);
}

static void check_inverse(MAT *A) {
  int n = A->m;
  MAT *Ainv = m_get(n, n), *I = m_get(n, n);
  m_inverse_small(A, Ainv);
  m_mlt(A, Ainv, I);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) CHECK_CLOSE(I->me[i][j], i == j ? 1.0 : 0.0, 1e-12);
  /* agrees with Meschach's general LU-based inverse */
  MAT *Aref = m_inverse(A, MNULL);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) CHECK_CLOSE(Ainv->me[i][j], Aref->me[i][j], 1e-12);
}

TEST_CASE(inverse_small_2x2) { check_inverse(mat2(4, 7, 2, 6)); }

TEST_CASE(inverse_small_3x3) {
  const double a[9] = {2, -3, 1, 2, 0, -1, 1, 4, 5};
  check_inverse(mat3(a));
  const double F[9] = {1.2, 0.3, 0.0, -0.1, 0.9, 0.05, 0.0, 0.02, 1.1};
  check_inverse(mat3(F));
}

/* The 3x3 helpers must agree with Meschach's general routines for any
 * input, including matrices with zero entries (axisymmetric tensors are
 * block-structured; an initial stretch tensor is the identity). */
static const double DENSE[9] = {1.5, -0.2, 0.7, 0.3, 2.0, -1.1, 0.0, 0.4, 0.9};
static const double PERMUT[9] = {0, 1, 0, 1, 0, 0, 0, 0, 1};   /* permutation */
static const double IDENT[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
static const double SHEAR[9] = {0, 0.25, 0, 0.25, 0, 0, 0, 0, 0}; /* like D */

static void check_mlt_3x3(const double a[9], const double b[9]) {
  MAT *A = mat3(a), *B = mat3(b), *C = m_get(3, 3), *Cref = m_get(3, 3);
  m_mlt_3x3(A, B, C);
  m_mlt(A, B, Cref);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) CHECK_CLOSE(C->me[i][j], Cref->me[i][j], TOL);
}

TEST_CASE(m_mlt_3x3_matches_m_mlt) {
  check_mlt_3x3(DENSE, DENSE);
  check_mlt_3x3(PERMUT, DENSE);
  check_mlt_3x3(DENSE, PERMUT);
  check_mlt_3x3(IDENT, SHEAR); /* V D with V = I at the first step */
  check_mlt_3x3(SHEAR, DENSE);
}

TEST_CASE(m_add_sub_3x3_match_meschach) {
  MAT *A = mat3(DENSE), *B = mat3(PERMUT), *C = m_get(3, 3), *R = m_get(3, 3);
  m_add_3x3(A, B, C);
  m_add(A, B, R);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) CHECK_CLOSE(C->me[i][j], R->me[i][j], TOL);
  m_sub_3x3(A, B, C);
  m_sub(A, B, R);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) CHECK_CLOSE(C->me[i][j], R->me[i][j], TOL);
}

TEST_CASE(contraction_is_frobenius_inner_product) {
  MAT *A = mat3(DENSE), *B = mat3(PERMUT);
  double ref = 0;
  for (int i = 0; i < 9; ++i) ref += DENSE[i] * PERMUT[i];
  CHECK_CLOSE(contraction(A, B), ref, TOL);
  /* non-square tensors (e.g. 2x3) */
  MAT *P = m_get(2, 3), *Q = m_get(2, 3);
  double r2 = 0;
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j) {
      P->me[i][j] = 1 + i + 2 * j;
      Q->me[i][j] = 0.5 * j - i;
      r2 += P->me[i][j] * Q->me[i][j];
    }
  CHECK_CLOSE(contraction(P, Q), r2, TOL);
}

TEST_CASE(smoothstep_endpoints_and_symmetry) {
  CHECK_CLOSE(smoothstep(-1.0, 2.0, 0.0), 0.0, TOL);
  CHECK_CLOSE(smoothstep(0.0, 2.0, 0.0), 0.0, TOL);
  CHECK_CLOSE(smoothstep(2.0, 2.0, 0.0), 1.0, TOL);
  CHECK_CLOSE(smoothstep(5.0, 2.0, 0.0), 1.0, TOL);
  CHECK_CLOSE(smoothstep(1.0, 2.0, 0.0), 0.5, TOL);
  CHECK_CLOSE(smoothstep(0.5, 2.0, 0.0) + smoothstep(1.5, 2.0, 0.0), 1.0, TOL);
}

/* ---------------------------------------------------------------- eig3 -- */

/* Checks A q_k = w_k q_k, Q orthonormal, and sum(w) = trace(A). */
static void check_eigensystem(const double a[9], int (*solver)(double **, double **, double *)) {
  double Abuf[3][3], Qbuf[3][3], w[3];
  double *A[3] = {Abuf[0], Abuf[1], Abuf[2]};
  double *Q[3] = {Qbuf[0], Qbuf[1], Qbuf[2]};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) Abuf[i][j] = a[3 * i + j];

  CHECK(solver(A, Q, w) == 0);

  for (int k = 0; k < 3; ++k)
    for (int i = 0; i < 3; ++i) {
      double Aq = 0;
      for (int j = 0; j < 3; ++j) Aq += a[3 * i + j] * Q[j][k];
      CHECK_CLOSE(Aq, w[k] * Q[i][k], 1e-9);
    }
  for (int k = 0; k < 3; ++k)
    for (int l = 0; l < 3; ++l) {
      double d = 0;
      for (int i = 0; i < 3; ++i) d += Q[i][k] * Q[i][l];
      CHECK_CLOSE(d, k == l ? 1.0 : 0.0, 1e-9);
    }
  CHECK_CLOSE(w[0] + w[1] + w[2], a[0] + a[4] + a[8], 1e-10);
}

static const double SYM_A[9] = {4, 1, -2, 1, 2, 0, -2, 0, 3};
static const double SYM_DIAG[9] = {5, 0, 0, 0, -1, 0, 0, 0, 2};
static const double SYM_REPEATED[9] = {2, 1, 1, 1, 2, 1, 1, 1, 2}; /* eigs 4,1,1 */
static const double SYM_STRETCH[9] = {1.44, 0.12, 0, 0.12, 0.81, 0, 0, 0, 1.0}; /* like C=F^T F */

TEST_CASE(dsyevh3_hybrid_solver) {
  check_eigensystem(SYM_A, dsyevh3);
  check_eigensystem(SYM_DIAG, dsyevh3);
  check_eigensystem(SYM_REPEATED, dsyevh3);
  check_eigensystem(SYM_STRETCH, dsyevh3);
}

TEST_CASE(dsyevq3_ql_solver) {
  check_eigensystem(SYM_A, dsyevq3);
  check_eigensystem(SYM_DIAG, dsyevq3);
  check_eigensystem(SYM_REPEATED, dsyevq3);
  check_eigensystem(SYM_STRETCH, dsyevq3);
}

TEST_CASE(dsyevc3_analytic_eigenvalues) {
  double Abuf[3][3], w[3];
  double *A[3] = {Abuf[0], Abuf[1], Abuf[2]};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) Abuf[i][j] = SYM_REPEATED[3 * i + j];
  dsyevc3(A, w);
  /* order is unspecified: compare sorted */
  for (int i = 0; i < 3; ++i)
    for (int j = i + 1; j < 3; ++j)
      if (w[j] < w[i]) { double t = w[i]; w[i] = w[j]; w[j] = t; }
  CHECK_CLOSE(w[0], 1.0, 1e-10);
  CHECK_CLOSE(w[1], 1.0, 1e-10);
  CHECK_CLOSE(w[2], 4.0, 1e-10);
}

int main(void) {
  RUN(determinant_2x2_and_3x3);
  RUN(trace_3x3);
  RUN(inverse_small_2x2);
  RUN(inverse_small_3x3);
  RUN(m_mlt_3x3_matches_m_mlt);
  RUN(m_add_sub_3x3_match_meschach);
  RUN(contraction_is_frobenius_inner_product);
  RUN(smoothstep_endpoints_and_symmetry);
  RUN(dsyevh3_hybrid_solver);
  RUN(dsyevq3_ql_solver);
  RUN(dsyevc3_analytic_eigenvalues);
  return TEST_SUMMARY();
}
