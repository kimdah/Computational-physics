# Simulation of a Penning trap (Project 3)
-------------------------------------------------
A Penning trap is a device that utilizes a static configuration of electic and magnetix fields to confine charged particles.
To build the code:
$ make

To run the code:
$ make run

To make plots:
$ make plot

![schematic penning trap](./figures/penning_trap.jpg)

# Folder Structure
Below, you will find a description of each folder. At the bottom you will find instructions on how to compile the program, run it and plot the datafiles that are produced.
## datafiles/
  This folder contains all the different datafiles(.txt files) that are produced when running the simulation. If the folder is empty, make sure you have compiled and run the program.

## figures/
  This folder will contain all the plots that are present in the report. If the folder is empty, make sure you ran the plot command in the terminal.

## include/
  This folder contains the header files for our program

## plot_files/
  This folder contains the python scrips that were used to plot the results. All the python files saves the figures in the 'figure/' folder. Note that if you try to run the program manually, it won't find the correct path to the datafiles, so make sure you run them using the '$make plot' command described below.

## src/
  This folder contains the Ising.cpp source file. Note that the main.cpp source file will be linked with the Ising.cpp source file.

The folder *figures* contains all figures from the .py files. \
The folder *datafiles* contains all output files from main.cpp.
