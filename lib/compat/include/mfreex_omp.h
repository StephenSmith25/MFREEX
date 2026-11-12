/*
 * mfreex_omp.h -- OpenMP compatibility shim.
 *
 * Include this instead of <omp.h>. When the compiler is building with
 * OpenMP (_OPENMP defined) it simply forwards to <omp.h>; otherwise it
 * provides serial stand-ins for the handful of runtime calls MFREEX uses,
 * so the code still builds and runs (single-threaded) on toolchains without
 * OpenMP, e.g. Apple Clang without libomp.
 */
#ifndef MFREEX_OMP_H_
#define MFREEX_OMP_H_

#ifdef _OPENMP
#include <omp.h>
#else
static inline void omp_set_num_threads(int n) { (void)n; }
static inline int omp_get_thread_num(void) { return 0; }
static inline int omp_get_num_threads(void) { return 1; }
static inline int omp_get_max_threads(void) { return 1; }
#endif

#endif /* MFREEX_OMP_H_ */
