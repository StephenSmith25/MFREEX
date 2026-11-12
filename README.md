<h1> MfreeX: A meshfree library written in C for non-linear explicit dynamics </h1>

This repositorya implementents a variant on the element-free Galerkin method. The current release of this library allows for the simulation of nonlinear 2D and axisymmetric problems using an explicit solver.

<h2> Features: </h2>

- Non linear material models (Hyperelastic, plastic, viscoelastic)
- Contact modelling.
- Nodal and Gauss integration.
- Clipped voronoi generator for implementing nodal integration.
- Implementation of both a total and updated element-free Galerkin method

<h2> To run the examples: </h2>

- Edit the root CMakeLists file to choose the problems to build 

<h2> Testing: </h2>

```bash
cmake -S . -B build && cmake --build build
ctest --test-dir build --output-on-failure      # all tests (~2 s)
ctest --test-dir build -L unit                  # unit tests only
```

- `unit_*`: property tests of the EFG building blocks. They cover MLS reproduction and derivative consistency, the SCNI/MSCNI patch tests, support domains versus brute-force neighbour search, the kernel functions, the 3x3 linear algebra and eigen-solvers, polar decomposition, and hyperelastic energy consistency (P = dW/dF).
- `golden_*`: rerun the `biaxial` (Buckley model) and `beam` (meshing) problems and compare their output with reference results in `tests/golden/`.
- `memcheck_*`: unit tests run under valgrind, if it is installed.

<h2> To do: </h2>

- Create input deck 
- Clean up redundant and commented out code.
- Continue user guide
- Lots of refactoring

