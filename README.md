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

<h2> To do: </h2>

- Create input deck 
- Clean up redundant and commented out code.
- Continue user guide
- Lots of refactoring

