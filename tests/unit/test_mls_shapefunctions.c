/*
 * Moving least squares (EFG) shape functions -- consistency tests.
 *
 * MLS shape functions built on a polynomial basis p(x) exactly reproduce
 * every polynomial in that basis:
 *      sum_I phi_I(x) p(x_I)        = p(x)
 *      sum_I dphi_I/dx_k(x) p(x_I)  = dp/dx_k(x)
 * for any node layout and weight function. These identities hold to
 * round-off, so any change to the shape-function code that breaks them is
 * a genuine regression.
 */
#include <stdlib.h>

#include "matrix.h"
#include "matrix2.h"
#include "mls_shapefunction.h"
#include "setDomain.h"

#include "mfreex_test.h"

#define TOL 1e-9

/* n x n grid on [0,1]^2, optionally perturbed so the layout is irregular */
static MAT *make_nodes(int n, double perturb) {
  MAT *nodes = m_get(n * n, 2);
  double h = 1.0 / (n - 1);
  srand(12345);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      double dx = 0, dy = 0;
      if (perturb > 0 && i > 0 && i < n - 1 && j > 0 && j < n - 1) {
        dx = perturb * h * (2.0 * rand() / RAND_MAX - 1.0);
        dy = perturb * h * (2.0 * rand() / RAND_MAX - 1.0);
      }
      nodes->me[i * n + j][0] = j * h + dx;
      nodes->me[i * n + j][1] = i * h + dy;
    }
  return nodes;
}

/* Evaluation points: interior, near edges and exactly at a corner */
static MAT *make_eval_points(void) {
  static const double pts[][2] = {{0.5, 0.5},   {0.123, 0.456}, {0.91, 0.07},
                                  {0.33, 0.99}, {0.0, 0.0},     {1.0, 0.5},
                                  {0.75, 0.25}, {0.051, 0.949}};
  int n = sizeof(pts) / sizeof(pts[0]);
  MAT *X = m_get(n, 2);
  for (int i = 0; i < n; ++i) {
    X->me[i][0] = pts[i][0];
    X->me[i][1] = pts[i][1];
  }
  return X;
}

/* Monomials of the basis, and their x/y derivatives */
static int n_monomials(const char *basis) {
  return basis[0] == 'q' ? 6 : 3;
}
static void monomials(const double *x, double *p, double *px, double *py) {
  double X = x[0], Y = x[1];
  double v[6] = {1, X, Y, X * X, X * Y, Y * Y};
  double dx[6] = {0, 1, 0, 2 * X, Y, 0};
  double dy[6] = {0, 0, 1, 0, X, 2 * Y};
  for (int k = 0; k < 6; ++k) {
    p[k] = v[k];
    px[k] = dx[k];
    py[k] = dy[k];
  }
}

static void check_reproduction(const char *basis, const char *weight,
                               double dmax, double perturb) {
  MAT *nodes = make_nodes(11, perturb);
  MAT *X = make_eval_points();
  double dmax_tensor[2] = {dmax, dmax};

  meshfreeDomain mfree = {.nodes = nodes,
                          .di = v_get(nodes->m),
                          .num_nodes = nodes->m,
                          .dim = 2,
                          .IS_AXI = 0,
                          .weight_function = (char *)weight,
                          .kernel_shape = "radial",
                          .basis_type = (char *)basis,
                          .is_constant_support_size = 1,
                          .dmax_radial = dmax,
                          .dmax_tensor = dmax_tensor,
                          .beta = 1.0};
  setDomain(&mfree);

  shape_function_container *sfc = mls_shapefunction(X, 2, &mfree);
  CHECK(sfc != NULL);
  CHECK(sfc->number_of_points == (int)X->m);

  int nm = n_monomials(basis);
  for (int i = 0; i < (int)X->m; ++i) {
    shape_function *sf = sfc->sf_list[i];
    double p[6], px[6], py[6];
    double sum[6] = {0}, sumx[6] = {0}, sumy[6] = {0};

    for (int k = 0; k < (int)sf->neighbours->max_dim; ++k) {
      double pI[6], dum1[6], dum2[6];
      monomials(nodes->me[sf->neighbours->ive[k]], pI, dum1, dum2);
      for (int m = 0; m < nm; ++m) {
        sum[m] += sf->phi->ve[k] * pI[m];
        sumx[m] += sf->dphi->me[k][0] * pI[m];
        sumy[m] += sf->dphi->me[k][1] * pI[m];
      }
    }

    monomials(X->me[i], p, px, py);
    for (int m = 0; m < nm; ++m) {
      CHECK_CLOSE(sum[m], p[m], TOL);
      CHECK_CLOSE(sumx[m], px[m], TOL);
      CHECK_CLOSE(sumy[m], py[m], TOL);
    }
  }
}

TEST_CASE(linear_basis_cubic_weight_regular_grid) {
  check_reproduction("linear", "cubic", 2.0, 0.0);
}
TEST_CASE(linear_basis_quartic_weight_regular_grid) {
  check_reproduction("linear", "quartic", 2.0, 0.0);
}
TEST_CASE(linear_basis_cubic_weight_irregular_nodes) {
  check_reproduction("linear", "cubic", 2.0, 0.3);
}
TEST_CASE(quadratic_basis_cubic_weight_regular_grid) {
  check_reproduction("quadratic", "cubic", 3.0, 0.0);
}
TEST_CASE(quadratic_basis_quartic_weight_irregular_nodes) {
  check_reproduction("quadratic", "quartic", 3.0, 0.3);
}

/* Shape functions evaluated at the nodes of a mesh should be smooth,
 * non-singular and (for EFG) not interpolating: phi_I(x_J) != delta_IJ in
 * general, but the sum is still 1. */
TEST_CASE(evaluation_at_nodes_partition_of_unity) {
  MAT *nodes = make_nodes(6, 0.2);
  double dmax_tensor[2] = {2.0, 2.0};
  meshfreeDomain mfree = {.nodes = nodes,
                          .di = v_get(nodes->m),
                          .num_nodes = nodes->m,
                          .dim = 2,
                          .weight_function = "cubic",
                          .kernel_shape = "radial",
                          .basis_type = "linear",
                          .is_constant_support_size = 1,
                          .dmax_radial = 2.0,
                          .dmax_tensor = dmax_tensor,
                          .beta = 1.0};
  setDomain(&mfree);
  shape_function_container *sfc = mls_shapefunction(nodes, 1, &mfree);
  for (int i = 0; i < (int)nodes->m; ++i) {
    shape_function *sf = sfc->sf_list[i];
    double s = 0;
    for (int k = 0; k < (int)sf->neighbours->max_dim; ++k) {
      CHECK(isfinite(sf->phi->ve[k]));
      s += sf->phi->ve[k];
    }
    CHECK_CLOSE(s, 1.0, TOL);
  }
}

int main(void) {
  RUN(linear_basis_cubic_weight_regular_grid);
  RUN(linear_basis_quartic_weight_regular_grid);
  RUN(linear_basis_cubic_weight_irregular_nodes);
  RUN(quadratic_basis_cubic_weight_regular_grid);
  RUN(quadratic_basis_quartic_weight_irregular_nodes);
  RUN(evaluation_at_nodes_partition_of_unity);
  return TEST_SUMMARY();
}
