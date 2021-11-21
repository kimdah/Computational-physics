Simulation of the two-dimensional Ising Model (Project 3)
-------------------------------------------------
The model will be used to explore temperature-dependent behaviour in ferromagnets and numerically estimate the critical temperature at which our 2D system undergoes a phase transition from a magnetized phase to a phase with no net magnetization.

The lattice of spins is defined by the state configuration matrix $s$.
...

We estimated the expectancy value of the energy per spin, $<eps> = $, and of the absolute value of the magnetization per spin, $<|m|>$. This were done by sampling $s$ using the Monte Carlo Markow Chain (MCMC) method and taking the mean of $epsilon$ and $|m|$ over MCMC cycles.
These were the used to find estimations for the specific heat capacity, $C__v$, and the susceptibility, $\chi$. These were both normalized to per spin.

We parallelized our code using OpenMP to decrease runtime for large state configuration matrices.  

-------------------
To build the code:  
$ make

To run the code:  
$ make run

To make plots:  
$ make plot

To test the code using chosen parameters:
$ ./main.exe <temperature (float)> <lattice side size (integer)> <MCMC cycles (integer)> <unordered lattice: use 0, ordered lattice: use -1 or 1> <output_file_name>
