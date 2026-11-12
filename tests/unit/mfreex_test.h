/*
 * mfreex_test.h -- minimal, dependency-free test helpers for MFREEX.
 *
 *   TEST_CASE(name) { ... CHECK(cond); CHECK_CLOSE(a, b, tol); ... }
 *   int main(void) { RUN(name); ...; return TEST_SUMMARY(); }
 */
#ifndef MFREEX_TEST_H_
#define MFREEX_TEST_H_

#include <math.h>
#include <stdio.h>

static int mfx_checks = 0;
static int mfx_failures = 0;

#define TEST_CASE(name) static void name(void)

#define RUN(name)                                                              \
  do {                                                                         \
    int before_ = mfx_failures;                                                \
    name();                                                                    \
    printf("[%s] %s\n", mfx_failures == before_ ? " OK " : "FAIL", #name);     \
  } while (0)

#define CHECK(cond)                                                            \
  do {                                                                         \
    ++mfx_checks;                                                              \
    if (!(cond)) {                                                             \
      ++mfx_failures;                                                          \
      fprintf(stderr, "  %s:%d: CHECK(%s) failed\n", __FILE__, __LINE__,       \
              #cond);                                                          \
    }                                                                          \
  } while (0)

/* |a - b| <= tol * max(1, |a|, |b|)  (relative for large values, absolute
 * near zero) */
#define CHECK_CLOSE(a, b, tol)                                                 \
  do {                                                                         \
    double a_ = (a), b_ = (b);                                                 \
    double s_ = fmax(1.0, fmax(fabs(a_), fabs(b_)));                           \
    ++mfx_checks;                                                              \
    if (!(fabs(a_ - b_) <= (tol) * s_)) {                                      \
      ++mfx_failures;                                                          \
      fprintf(stderr, "  %s:%d: CHECK_CLOSE(%s, %s) failed: %.17g vs %.17g\n", \
              __FILE__, __LINE__, #a, #b, a_, b_);                             \
    }                                                                          \
  } while (0)

#define TEST_SUMMARY()                                                         \
  (printf("%d checks, %d failures\n", mfx_checks, mfx_failures),              \
   mfx_failures == 0 ? 0 : 1)

#endif /* MFREEX_TEST_H_ */
