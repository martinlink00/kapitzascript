# Kapitza Script
The scripts in this directory use the method of modulated Fourier expansion to calculate the effective potential governing the slow oscillation mode of an inverted or Kapitza pendulum in an iterative calculation using Mathematica 12.1.
# How to use the scripts
  - Calculate potential corrections by typing math -script potential.wls in the project folder
  - Follow command line instructions to input parameters
  - After the potential script has successfully run, match initial conditions by running math -script incond.wls
  - Results as well as logfiles will be saved into a subdirectory exports, which will be automatically created when the program is first run
# Results
  - View results using the Results.nb notebook in the project folder (instructions on usage given there)
  - Each new run of either one of the scripts will overwrite previous results!
    - If you wish to keep results, copy the exports directory and Results.nb somewhere else to avoid risk of loosing everything 
