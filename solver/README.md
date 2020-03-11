# Finite volume method for shallow water equations

The code solves the shallow water equations for a dam break problem using the finite volume method. The the Lax–Friedrichs or Roe methods are used to evaluate the fluxes. The Euler method is used for the explicit time stepping.
The code is for teaching/explanation purposes and is not optimised in any way. 

# Compile

The code requires the Eigen C++ libary. Download the libary and extract its contents to a folder named Eigen inside the includes folder. Once done run make to compile the code. The sw_solver binary will be created in the bin folder. 

# Run

To run the code navigate to the folder of the case e.g. examples/dam_lf and run ./sw_solver dam. The program will solve the dam break problem using the Lax-Friedrichs method for 500 time steps.

The sw_solver binary requires the following files:

- .conf - contains the solver parameters such as timestep, number of timesteps, Riemann solver, etc.
- .grid - contains the grid