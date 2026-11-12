/*
 * EFG nodal integration (SCNI / MSCNI) -- patch tests.
 *
 * Stabilised conforming nodal integration replaces the strain at a node by
 * its average over the node's (clipped) Voronoi cell, computed via the
 * divergence theorem as a boundary integral of the MLS shape functions.
 * Because MLS with a linear basis reproduces linear fields exactly, SCNI
 * must reproduce a homogeneous deformation u(x) = A x + c exactly in every
 * cell (the linear patch test), and the cells must tile the domain.
 *
 * B-matrix rows (generate_Bmat, 2D): [du/dx, dv/dy, du/dy, dv/dx(, u/r)]
 * acting on the interleaved displacement vector [u_I, v_I, ...].
 */
#include <stdlib.h>

#include "matrix.h"
#include "matrix2.h"

#include "Integration/SCNI/generate_mscni.h"
#include "Integration/SCNI/generate_scni.h"
#include "generate_voronoi.h"
#include "mls_shapefunction.h"
#include "setDomain.h"

#include "mfreex_test.h"

#define TOL 1e-9

typedef struct {
  MAT *nodes;
  int *boundary;
  int num_boundary;
  meshfreeDomain mfree;
  double dmax_tensor[2];
  voronoi_diagram *vor;
} patch;

/* n x n nodes on [x0, x0+L] x [0, L]; interior nodes optionally perturbed.
 * Boundary node indices listed counter-clockwise. */
static patch *make_patch(int n, double x0, double L, double perturb) {
  patch *p = calloc(1, sizeof(patch));
  double h = L / (n - 1);
  p->nodes = m_get(n * n, 2);
  srand(4242);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      double dx = 0, dy = 0;
      if (perturb > 0 && i > 0 && i < n - 1 && j > 0 && j < n - 1) {
        dx = perturb * h * (2.0 * rand() / RAND_MAX - 1.0);
        dy = perturb * h * (2.0 * rand() / RAND_MAX - 1.0);
      }
      p->nodes->me[i * n + j][0] = x0 + j * h + dx;
      p->nodes->me[i * n + j][1] = i * h + dy;
    }

  p->boundary = malloc(4 * (n - 1) * sizeof(int));
  int b = 0;
  for (int j = 0; j < n - 1; ++j) p->boundary[b++] = j;                      /* bottom */
  for (int i = 0; i < n - 1; ++i) p->boundary[b++] = i * n + (n - 1);        /* right  */
  for (int j = n - 1; j > 0; --j) p->boundary[b++] = (n - 1) * n + j;        /* top    */
  for (int i = n - 1; i > 0; --i) p->boundary[b++] = i * n;                  /* left   */
  p->num_boundary = b;

  p->dmax_tensor[0] = p->dmax_tensor[1] = 2.0;
  meshfreeDomain m = {.nodes = p->nodes,
                      .di = v_get(p->nodes->m),
                      .num_nodes = p->nodes->m,
                      .dim = 2,
                      .weight_function = "cubic",
                      .kernel_shape = "radial",
                      .basis_type = "linear",
                      .is_constant_support_size = 1,
                      .dmax_radial = 2.0,
                      .dmax_tensor = p->dmax_tensor,
                      .beta = 1.0};
  p->mfree = m;
  setDomain(&p->mfree);
  p->vor = generate_voronoi(p->nodes->base, p->boundary, p->nodes->m, p->num_boundary, 2);
  return p;
}

/* u = A x + c */
static const double A[2][2] = {{0.013, -0.007}, {0.021, 0.004}};
static const double C[2] = {0.3, -0.1};

static double cell_strain(MAT *B, IVEC *idx, MAT *nodes, int row) {
  double s = 0;
  for (int k = 0; k < (int)idx->max_dim; ++k) {
    double *x = nodes->me[idx->ive[k]];
    double u = A[0][0] * x[0] + A[0][1] * x[1] + C[0];
    double v = A[1][0] * x[0] + A[1][1] * x[1] + C[1];
    s += B->me[row][2 * k] * u + B->me[row][2 * k + 1] * v;
  }
  return s;
}

/* Clipped Voronoi cells must be simple polygons: no vertex repeated
 * consecutively (a zero-length edge gives 0/0 in the SCNI normal, and skews
 * the vertex-average "centre" used as r in axisymmetric problems). */
static void check_voronoi_cells(patch *p) {
  voronoi_diagram *v = p->vor;
  CHECK(v->num_cells == (int)p->nodes->m);
  for (int i = 0; i < v->num_cells; ++i) {
    int m = v->num_cell_verticies[i];
    CHECK(m >= 3);
    for (int k = 0; k < m; ++k)
      CHECK(v->index[i][k] != v->index[i][(k + 1) % m]);
  }
}

TEST_CASE(voronoi_cells_have_no_repeated_vertices) {
  check_voronoi_cells(make_patch(9, 0.0, 1.0, 0.0));
  check_voronoi_cells(make_patch(9, 0.0, 1.0, 0.25));
  check_voronoi_cells(make_patch(7, 0.5, 1.0, 0.0));
}

static void check_scni_patch(double perturb) {
  patch *p = make_patch(9, 0.0, 1.0, perturb);
  SCNI_OBJ *s = generate_scni(p->vor, NULL, 0, 0, 2, &p->mfree);
  CHECK(s->num_points == (int)p->nodes->m);

  double area = 0;
  for (int i = 0; i < s->num_points; ++i) {
    SCNI *c = s->scni[i];
    area += c->area;
    CHECK_CLOSE(cell_strain(c->B, c->sfIndex, p->nodes, 0), A[0][0], TOL); /* du/dx */
    CHECK_CLOSE(cell_strain(c->B, c->sfIndex, p->nodes, 1), A[1][1], TOL); /* dv/dy */
    CHECK_CLOSE(cell_strain(c->B, c->sfIndex, p->nodes, 2), A[0][1], TOL); /* du/dy */
    CHECK_CLOSE(cell_strain(c->B, c->sfIndex, p->nodes, 3), A[1][0], TOL); /* dv/dx */
  }
  CHECK_CLOSE(area, 1.0, 1e-12); /* Voronoi cells tile the unit square */
}

TEST_CASE(scni_linear_patch_test_regular) { check_scni_patch(0.0); }
TEST_CASE(scni_linear_patch_test_irregular) { check_scni_patch(0.25); }

/* Axisymmetric: the extra B row is the hoop strain u_r / r. For
 * u_r = a r + b z the hoop strain at node I is (a r_I + b z_I) / r_I. */
TEST_CASE(scni_axisymmetric_hoop_strain) {
  patch *p = make_patch(7, 0.5, 1.0, 0.0); /* r in [0.5, 1.5]: no r = 0 */
  SCNI_OBJ *s = generate_scni(p->vor, NULL, 0, 1, 2, &p->mfree);
  for (int i = 0; i < s->num_points; ++i) {
    SCNI *c = s->scni[i];
    CHECK(c->B->m == 5);
    double hoop = 0;
    for (int k = 0; k < (int)c->sfIndex->max_dim; ++k) {
      double *x = p->nodes->me[c->sfIndex->ive[k]];
      double ur = A[0][0] * x[0] + A[0][1] * x[1] + C[0];
      hoop += c->B->me[4][2 * k] * ur;
    }
    double *xi = p->nodes->me[i];
    CHECK_CLOSE(hoop, (A[0][0] * xi[0] + A[0][1] * xi[1] + C[0]) / xi[0], TOL);
    /* in-plane rows must be unaffected by the axisymmetric option */
    CHECK_CLOSE(cell_strain(c->B, c->sfIndex, p->nodes, 0), A[0][0], TOL);
  }
}

TEST_CASE(mscni_linear_patch_test) {
  patch *p = make_patch(9, 0.0, 1.0, 0.25);
  MSCNI_OBJ *s = generate_mscni(p->vor, NULL, 0, 0, &p->mfree);
  CHECK(s->num_points == (int)p->nodes->m);
  double area = 0;
  for (int i = 0; i < s->num_points; ++i) {
    MSCNI *c = s->scni[i];
    area += c->area;
    CHECK_CLOSE(cell_strain(c->B, c->sfIndex, p->nodes, 0), A[0][0], TOL);
    CHECK_CLOSE(cell_strain(c->B, c->sfIndex, p->nodes, 1), A[1][1], TOL);
    CHECK_CLOSE(cell_strain(c->B, c->sfIndex, p->nodes, 2), A[0][1], TOL);
    CHECK_CLOSE(cell_strain(c->B, c->sfIndex, p->nodes, 3), A[1][0], TOL);
  }
  CHECK_CLOSE(area, 1.0, 1e-12);
}

int main(void) {
  RUN(voronoi_cells_have_no_repeated_vertices);
  RUN(scni_linear_patch_test_regular);
  RUN(scni_linear_patch_test_irregular);
  RUN(scni_axisymmetric_hoop_strain);
  RUN(mscni_linear_patch_test);
  return TEST_SUMMARY();
}
