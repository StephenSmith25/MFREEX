/*
 * Known bugs in master, recorded as tests.
 *
 * Each case is registered with CTest as WILL_FAIL, so the suite stays green
 * while the bug exists. When a bug is fixed, CTest reports the test as
 * failing ("passed unexpectedly") -- remove WILL_FAIL for that case in
 * tests/CMakeLists.txt so it becomes an ordinary regression test.
 *
 * Usage: test_known_bugs <case>
 */
#include <string.h>

#include "matrix.h"
#include "matrix2.h"

#include "m_inverse_small.h"

#include "mfreex_test.h"

int SVK(VEC *stressVoigt, MAT *defGrad, VEC *params);

/* SVK in 3D builds the Green-Lagrange strain from F F^T (rows) instead of
 * F^T F (columns). For a plane deformation (F13=F23=F31=F32=0, F33=1) the
 * in-plane 1st PK stress must equal the 2D result; under simple shear it
 * does not. */
static int svk_3d_matches_2d_under_shear(void) {
  double g = 0.3;
  VEC *params = v_get(2);
  params->ve[0] = 1200.0;
  params->ve[1] = 800.0;

  MAT *F2 = m_get(2, 2);
  F2->me[0][0] = 1; F2->me[0][1] = g; F2->me[1][1] = 1;
  VEC *P2 = v_get(4);
  SVK(P2, F2, params);

  MAT *F3 = m_get(3, 3);
  F3->me[0][0] = 1; F3->me[0][1] = g; F3->me[1][1] = 1; F3->me[2][2] = 1;
  VEC *P3 = v_get(9);
  SVK(P3, F3, params);

  /* 3D Voigt ordering from SVK.c: [P11, P22, P33, ...]; compare P11, P22 */
  CHECK_CLOSE(P3->ve[0], P2->ve[0], 1e-12);
  CHECK_CLOSE(P3->ve[1], P2->ve[1], 1e-12);
  return TEST_SUMMARY();
}

/* m_inverse_small(A, NULL) is meant to allocate the output, but the
 * allocation is assigned to a shadowing local, so it dereferences NULL. */
static int m_inverse_small_allocates_output(void) {
  MAT *A = m_get(2, 2);
  A->me[0][0] = 4; A->me[0][1] = 7; A->me[1][0] = 2; A->me[1][1] = 6;
  MAT *Ainv = m_inverse_small(A, MNULL);
  CHECK(Ainv != MNULL);
  return TEST_SUMMARY();
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s <case>\n", argv[0]);
    return 2;
  }
  if (strcmp(argv[1], "svk_3d_shear") == 0) return svk_3d_matches_2d_under_shear();
  if (strcmp(argv[1], "m_inverse_small_null_out") == 0) return m_inverse_small_allocates_output();
  fprintf(stderr, "unknown case %s\n", argv[1]);
  return 2;
}
