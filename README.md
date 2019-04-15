<h1> MfreeX: A meshfree library written in C for non-linear explicit dynamics </h1>

This repository contains the code for a that implementents a variant on the element-free Galerkin method. The current release of this library allows for the simulation of nonlinear 2D and axisymmetric problems using an explicit solver.

<h2> Features: </h2>

- Non linear material models (Hyperelastic, plastic, viscoelastic)
- Contact modelling.
- Clipped voronoi generator for implementing nodal integration.

<h2> To run the examples: </h2>

- Edit the root CMakeLists file to choose the problems to build 

<h2> To do: </h2>

- Create input deck 
- Clean up redundant code 
- Continue user guide
- <h1> Lots of refactoring </h1>

<h1> Notes </h1>
This is currently a work in progress and therefore a number of issues are expected to exist within the code. 
