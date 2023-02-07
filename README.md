# FluidSim

To do:
 - Parallel computing.
 - - CUDA: to do.
 - Remove unused parts of Eigen library.
 - Correct handling of solid cells and domain boundaries.
 - Improve particle seeding. Intialise some grid cells as fluid, then seed particles. Rather than current method of placing particles.
 - Improve interpolation from grid to particles, B-spline as well?
 

Long term:
 - Application refactoring
 - - Increase modularity
 - Use templates alongside inheritance for 2D and 3D class variants.
 - Improve method of changing between 2D and 3D simulations.
