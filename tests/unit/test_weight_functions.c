/*
 * EFG kernel (weight) functions: cubic and quartic splines.
 *
 * w(r) on normalised distance r = |x - x_I| / d_I must be
 *   - positive inside the support, zero outside (compact support),
 *   - C1 continuous (C2 for the quartic) including at r = 1/2 and r = 1,
 *   - consistent: the returned derivative matches the function.
 * Weight arrays are [w, dw/dr, (d2w/dr2)].
 */
#include "matrix.h"
#include "matrix2.h"

#include "cubic_spline.h"
#include "quartic_spline.h"

#include "mfreex_test.h"

typedef int (*spline_fn)(VEC *, double);

static double W(spline_fn f, double r, int k) {
  VEC *w = v_get(3);
  f(w, r);
  double out = w->ve[k];
  V_FREE(w);
  return out;
}

static void check_spline(spline_fn f, int order) {
  /* compact support */
  CHECK(W(f, 0.0, 0) > 0);
  CHECK(W(f, 0.99, 0) > 0);
  CHECK_CLOSE(W(f, 1.0, 0), 0.0, 1e-12);
  CHECK_CLOSE(W(f, 1.5, 0), 0.0, 0.0);
  CHECK_CLOSE(W(f, 1.5, 1), 0.0, 0.0);
  /* zero slope at the centre and at the edge of the support */
  CHECK_CLOSE(W(f, 0.0, 1), 0.0, 1e-12);
  CHECK_CLOSE(W(f, 1.0, 1), 0.0, 1e-12);
  /* monotonically decreasing on (0, 1) */
  for (double r = 0.01; r < 1.0; r += 0.01) CHECK(W(f, r, 0) < W(f, r - 0.01, 0));
  /* derivative(s) consistent with the function (central differences) */
  const double h = 1e-6;
  for (double r = 0.05; r < 1.0; r += 0.0625) {
    double fd1 = (W(f, r + h, 0) - W(f, r - h, 0)) / (2 * h);
    CHECK_CLOSE(W(f, r, 1), fd1, 1e-6);
    if (order >= 2) {
      double fd2 = (W(f, r + h, 1) - W(f, r - h, 1)) / (2 * h);
      CHECK_CLOSE(W(f, r, 2), fd2, 1e-5);
    }
  }
  /* continuity across the piecewise break point r = 1/2 */
  const double e = 1e-10;
  CHECK_CLOSE(W(f, 0.5 - e, 0), W(f, 0.5 + e, 0), 1e-8);
  CHECK_CLOSE(W(f, 0.5 - e, 1), W(f, 0.5 + e, 1), 1e-8);
}

TEST_CASE(cubic_spline_properties) { check_spline(cubic_spline, 1); }
TEST_CASE(quartic_spline_properties) { check_spline(quartic_spline, 2); }

int main(void) {
  RUN(cubic_spline_properties);
  RUN(quartic_spline_properties);
  return TEST_SUMMARY();
}
