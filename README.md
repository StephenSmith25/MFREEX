# MFREEX

**A meshfree library written in C for non-linear explicit dynamics.**

MFREEX implements a variant of the element-free Galerkin (EFG) method for
non-linear 2D and axisymmetric problems using an explicit solver. It was
developed as part of a PhD project at Queen's University Belfast.

## Features

- Non-linear material models (hyperelastic, plastic, viscoelastic / Buckley)
- Contact modelling
- Nodal (SCNI / MSCNI) and Gauss integration
- Clipped Voronoi generator for nodal integration
- Total and updated Lagrangian EFG formulations

## Building

Requirements: a C99 compiler (GCC, Clang or Apple Clang) and CMake ≥ 3.16.
OpenMP is optional: without it MFREEX builds and runs single-threaded.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

With CMake ≥ 3.21 you can use the presets instead:

```bash
cmake --preset release && cmake --build --preset release   # optimised, OpenMP if found
cmake --preset debug   && cmake --build --preset debug     # -O0 -g, ASan + UBSan, serial
cmake --preset serial  && cmake --build --preset serial    # optimised, no OpenMP
```

### macOS

Apple Clang does not ship OpenMP. Either build serially (it just works), or
install `libomp` from Homebrew for a threaded build:

```bash
brew install libomp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DOpenMP_ROOT="$(brew --prefix libomp)"
cmake --build build --parallel
# or: cmake --preset macos && cmake --build --preset macos
```

### Build options

| Option | Default | Description |
| --- | --- | --- |
| `MFREEX_USE_OPENMP` | `ON` | Use OpenMP threading if available |
| `MFREEX_BUILD_PROBLEMS` | `ON` | Build the example problems |
| `MFREEX_ENABLE_WARNINGS` | `OFF` | `-Wall -Wextra` for the MFREEX library sources |
| `MFREEX_WARNINGS_AS_ERRORS` | `OFF` | `-Werror` for the MFREEX library sources |
| `MFREEX_SANITIZE` | `OFF` | AddressSanitizer + UndefinedBehaviorSanitizer |
| `PROBLEM_<NAME>` | see `CMakeLists.txt` | Enable/disable individual problems |

A `compile_commands.json` is generated in the build directory for clangd and
other editor tooling.

## Running the examples

Each problem is built into `build/bin/<problem>/` together with its input
files (`*.nodes`, `*.segs`, ...). The executables read inputs and write output
relative to the working directory, so run them from there:

```bash
cd build/bin/rubber
mkdir -p Displacement Stress
OMP_NUM_THREADS=4 ./rubber
```

Thread count is controlled with `OMP_NUM_THREADS`, though some routines still
set it explicitly (see *Known issues*).

Post-processing scripts are in `matlab/` and `python/`
(`pip install -r python/requirements.txt`).

## Repository layout

```
MFREEX/          core library (include/, src/)
problems/        example problems, one executable each
lib/             vendored dependencies (see below)
cmake/           CMake helper modules
input/           experimental input-deck parser
matlab/, python/ post-processing scripts and reference data
userguide/       LaTeX user guide (userguide.pdf)
```

To add a new problem, create `problems/<name>/` with your source and input
files, a one-line `CMakeLists.txt`:

```cmake
mfreex_add_problem(<name> SOURCES <name>.c)
```

and an `add_subdirectory()` entry in `problems/CMakeLists.txt`.

## Vendored dependencies and licences

All dependencies are bundled in `lib/` and built automatically; nothing needs
to be installed. Note that **some carry non-commercial licence terms**.

| Directory | CMake target | Library | Licence |
| --- | --- | --- | --- |
| `lib/Meschach` | `meschach` | Meschach matrix library (Stewart & Leyk) | Meschach licence: free to copy/modify; no charge may be made; changes must be noted |
| `lib/triangulate` | `triangulate` | Triangle 1.6 (J. R. Shewchuk) + wrappers | **Non-commercial use only** without permission |
| `lib/gpc` | `gpc` | General Polygon Clipper 2.32 (A. Murta) | **Non-commercial use only** without permission |
| `lib/Eigen` | `eig3` | 3×3 eigen-solvers (J. Kopp). *Not* the C++ Eigen library | LGPL-2.1+ |
| `lib/clipped_voronoi` | `CVLIB` | Clipped Voronoi, using `jc_voronoi` | MIT (jc_voronoi) |
| `lib/thread_pool` | `thread_pool` | C-Thread-Pool (Pithikos) | MIT |
| `lib/efsflib`, `lib/util` | `EFSFLIB`, `util` | EFG shape functions and helpers (part of this project) | — |

## Known issues

- Several example problems crash at runtime, in this version and in the
  original code. Running the `debug` preset (ASan/UBSan) shows where. Remaining
  reports include a stack overflow inside Triangle (`cylinder`), misaligned
  stores in `jc_voronoi` (`taylor_bar`, `rubber_comp`), a zero-length VLA in
  `rubber`, and a `memcpy` from `NULL` in `material_point.c`
  (`beam_updated`, `rigid_punch`).
- `problems/beam` calls `exit(0)` straight after meshing (left-over debug code).
- Some library routines hard-code `omp_set_num_threads(3|4)`, which overrides
  `OMP_NUM_THREADS`.

## To do

- Create input deck
- Clean up redundant and commented-out code
- Continue user guide
- Lots of refactoring
