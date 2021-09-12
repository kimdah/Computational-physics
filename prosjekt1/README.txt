Programs for Project 1: numerical solution to the one-dimensional Poisson equation by converting the second derivative to a solvable matrix equation Av = g.

main.cpp
--------
C++ code that computes the exact solution, u(x), and the numerical approximations v(x) and v^(x) using both the general algorithm and special algorithm 
for solving the matrix equation. This is done with different number of steps, N = 10**i for i = 1,2,...
Outputs the result to files based on N:
exact_data<N>.txt
approx_general<N>.txt
approx_special<N>.txt

Build commands: 
g++ -c -std=c++11 main.cpp
g++ main.o -o main.exe -larmadillo

Run command: ./main.exe N


plot_exact.py
-------
Python script that reads the data in ....

output.txt and generates plots of the absolute error, 
relative error and log10(relative error). Plots are saved as pdf files.

Run command: python3 plot.py
