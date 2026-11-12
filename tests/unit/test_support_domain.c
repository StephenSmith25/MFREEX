/*
 * EFG support domains (setDomain) and neighbour search (get_point_neighbours).
 *
 *  - Support sizes must be deterministic: they may not depend on stack
 *    contents left by earlier calls.
 *  - The method is unit-independent: scaling all coordinates by s scales
 *    the supports by s and leaves neighbour sets unchanged (models in mm
 *    vs m must behave identically).
 *  - Neighbour search must agree with a brute-force definition of the
 *    support (|x - x_I| < d_I radial; |x_k - x_Ik| <= d_Ik rectangular).
 */
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "matrix2.h"
#include "neighbours.h"
#include "setDomain.h"

#include "mfreex_test.h"

static MAT *grid(int n, double scale, double perturb) {
  MAT *X = m_get(n * n, 2);
  double h = 1.0 / (n - 1);
  srand(99);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      double dx = perturb * h * (2.0 * rand() / RAND_MAX - 1.0);
      double dy = perturb * h * (2.0 * rand() / RAND_MAX - 1.0);
      X->me[i * n + j][0] = scale * (j * h + dx);
      X->me[i * n + j][1] = scale * (i * h + dy);
    }
  return X;
}

static meshfreeDomain domain(MAT *X, const char *shape, double *dmax_tensor) {
  meshfreeDomain m = {.nodes = X,
                      .di = v_get(X->m),
                      .num_nodes = X->m,
                      .dim = 2,
                      .weight_function = "cubic",
                      .kernel_shape = (char *)shape,
                      .basis_type = "linear",
                      .is_constant_support_size = 0,
                      .dmax_radial = 2.0,
                      .dmax_tensor = dmax_tensor,
                      .beta = 1.0};
  return m;
}

/* Leave large garbage values in the stack region setDomain will reuse */
static void __attribute__((noinline)) dirty_stack(void) {
  volatile double junk[4096];
  for (int i = 0; i < 4096; ++i) junk[i] = 1e300;
}

TEST_CASE(rectangular_supports_are_deterministic) {
  MAT *X = grid(8, 1.0, 0.2);
  double dt[2] = {1.5, 1.5};

  meshfreeDomain a = domain(X, "rectangular", dt);
  setDomain(&a);
  dirty_stack();
  meshfreeDomain b = domain(X, "rectangular", dt);
  setDomain(&b);

  for (int i = 0; i < (int)X->m; ++i)
    for (int k = 0; k < 2; ++k) {
      CHECK(isfinite(b.di_tensor->me[i][k]) && b.di_tensor->me[i][k] < 10.0);
      CHECK_CLOSE(a.di_tensor->me[i][k], b.di_tensor->me[i][k], 0.0);
    }
}

static int same_set(IVEC *a, IVEC *b) {
  if (a->max_dim != b->max_dim) return 0;
  for (int i = 0; i < (int)a->max_dim; ++i)
    if (a->ive[i] != b->ive[i]) return 0;
  return 1;
}

static void check_scale_invariance(const char *shape, double scale) {
  double dt[2] = {1.5, 1.5};
  MAT *X1 = grid(8, 1.0, 0.2), *Xs = grid(8, scale, 0.2);
  meshfreeDomain m1 = domain(X1, shape, dt), ms = domain(Xs, shape, dt);
  setDomain(&m1);
  setDomain(&ms);
  static const double q[][2] = {{0.5, 0.5}, {0.137, 0.861}, {0.99, 0.02}, {0.3, 0.3}};
  for (int p = 0; p < 4; ++p) {
    double x1[2] = {q[p][0], q[p][1]}, xs[2] = {scale * q[p][0], scale * q[p][1]};
    IVEC *n1 = get_point_neighbours(iv_get(10), x1, &m1);
    IVEC *ns = get_point_neighbours(iv_get(10), xs, &ms);
    CHECK(n1->max_dim > 0);
    CHECK(same_set(n1, ns));
  }
}

TEST_CASE(radial_neighbours_independent_of_units) {
  check_scale_invariance("radial", 1000.0);
  check_scale_invariance("radial", 1e-3);
}
TEST_CASE(rectangular_neighbours_independent_of_units) {
  check_scale_invariance("rectangular", 1000.0);
  check_scale_invariance("rectangular", 1e-3);
}

static void check_brute_force(const char *shape) {
  double dt[2] = {1.5, 1.5};
  MAT *X = grid(8, 1.0, 0.2);
  meshfreeDomain m = domain(X, shape, dt);
  setDomain(&m);
  int rect = strcmp(shape, "rectangular") == 0;
  for (double qx = 0.03; qx < 1.0; qx += 0.11)
    for (double qy = 0.05; qy < 1.0; qy += 0.13) {
      double x[2] = {qx, qy};
      IVEC *nb = get_point_neighbours(iv_get(10), x, &m);
      int count = 0;
      for (int i = 0; i < (int)X->m; ++i) {
        double dx = fabs(x[0] - X->me[i][0]), dy = fabs(x[1] - X->me[i][1]);
        int inside = rect ? (dx <= m.di_tensor->me[i][0] && dy <= m.di_tensor->me[i][1])
                          : (sqrt(dx * dx + dy * dy) < m.di->ve[i]);
        if (inside) {
          CHECK(count < (int)nb->max_dim && nb->ive[count] == i);
          ++count;
        }
      }
      CHECK(count == (int)nb->max_dim);
    }
}

TEST_CASE(radial_neighbours_match_brute_force) { check_brute_force("radial"); }
TEST_CASE(rectangular_neighbours_match_brute_force) { check_brute_force("rectangular"); }

int main(void) {
  RUN(rectangular_supports_are_deterministic);
  RUN(radial_neighbours_independent_of_units);
  RUN(rectangular_neighbours_independent_of_units);
  RUN(radial_neighbours_match_brute_force);
  RUN(rectangular_neighbours_match_brute_force);
  return TEST_SUMMARY();
}
