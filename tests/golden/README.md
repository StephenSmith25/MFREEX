# Golden reference output

Reference results for the example problems that run to completion,
generated from `master` at f4264ba6 (the original thesis code) with
`OMP_NUM_THREADS=1`. The `golden_*` tests rerun each problem and compare
its output with `tests/tools/numdiff` (relative tolerance 1e-6).

| Directory | Problem | What it exercises |
| --- | --- | --- |
| `biaxial/` | `problems/Biaxial` | Buckley viscoelastic model (PET) under equal-biaxial stretch, 1875 load steps |
| `beam/` | `problems/beam` (`beam_triangle.c`) | Triangle meshing, boundary detection, background search cells |

If a deliberate change alters these results, regenerate the files by
running the problem and copying its output here, and say why in the commit
message.
