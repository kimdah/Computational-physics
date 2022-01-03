# Numerical solution to the one-dimensional Poisson equation by converting the second derivative to a solvable matrix equation Av = g (Project 1)
---------------------------

main.cpp
--------
C++ code that computes the exact solution, u(x), and the numerical approximations v(x) and v^(x) using both the general algorithm and special algorithm for solving the matrix equation.
This is done with different number of steps, N = 10**i for i = 1,2,...,6

The code outputs the result to files based on N:
exact_data<N>.txt
approx_general<N>.txt
approx_special<N>.txt

Stores log10 absolute and relative error, as well as max(rel_error) in file errors<N>.txt

plot_exact.py
-------
Python script that reads the data in exact_data{N}.txt and plots the exact solution of the Poisson equation, u(x).
Output: exact.pdf


general_vs_exact.py
-------
Plots the general approximations from approx_general<N>.txt against the exact values exact_data<N>.txt.
Output: general_vs_exact.pdf



errorplot.py
-------
Plots logarithmic10 absolute and relative errors from file errors<N>.txt against x-values.
Output:'error_plot_log.pdf'


# Compiling, running and plotting

To build the code:  
$ make

To run the code:  
$ make run

To make plots:  
$ make plot


# Folder Structure
Below, you will find a description of each folder. At the bottom you will find instructions on how to compile the program, run it and plot the datafiles that are produced.
## datafiles/
  This folder contains all the different datafiles(.txt files) that are produced when running the simulation. If the folder is empty, make sure you have compiled and run the program.

## figures/
  This folder will contain all the plots that are present in the report. If the folder is empty, make sure you ran the plot command in the terminal.

## plot_files/
  This folder contains the python scrips that were used to plot the results. All the python files saves the figures in the 'figure/' folder. Note that if you try to run the program manually, it won't find the correct path to the datafiles, so make sure you run them using the '$make plot' command described below.
