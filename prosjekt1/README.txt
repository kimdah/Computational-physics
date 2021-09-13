Programs for Project 1: numerical solution to the one-dimensional Poisson equation by converting the second derivative to a solvable matrix equation Av = g.

main.cpp
--------
C++ code that computes the exact solution, u(x), and the numerical approximations v(x) and v^(x) using both the general algorithm and special algorithm 
for solving the matrix equation. This is done with different number of steps, N = 10**i for i = 1,2,...
Outputs the result to files based on N:
exact_data<N>.txt
approx_general<N>.txt
approx_special<N>.txt

Stores log10 absolute and relative error, as well as max(rel_error) in fil errors<N>.txt

Build commands: 
g++ -c -std=c++11 main.cpp
g++ main.o -o main.exe -larmadillo

Run command: ./main.exe N


plot_exact.py
-------
Python script that reads the data in exact_data1000.txt and plots the exact solution of the Poisson equation, u(x). 
Output: problem2.pdf

Run command: python plot_exact.py

general_vs_exact.py
-------
Plots the general approximations from approx_general<N>.txt against the exact values exact_data<N>.txt.
Output: general_vs_exact.pdf

Run command: python general_vs_exact.py


errorplot.py
-------
Plots logarithmic10 absolute and relative errors from file errors<N>.txt against x-values.
Output:'error_plot_log.pdf'

Run command: python errorplot.py
