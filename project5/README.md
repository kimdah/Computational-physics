# Simulating the two-dimensional time-dependent Schrödinger equation
## Desctiption
This folder is contains source code, datafiles and figures for project 5 in FYS4150/3150. This project aims to study the behavoiur of a quantum particle in a box by using the normalizes schrodiner equation. We conduct several experiments using this simulation, for example the single-, double- and triple-slit experiment. Instuctions on compiling and running the function, as well as plotting the figures can be found at the bottom of this README file. 

## Folder structure

### `src/`
This folder contains the source code `Crank.cpp` to our program, which is linked during compilation. 

### `plotfiles/`
This folder contains the python scrips that are used to plot the animations and figures we have used in the report. Some of these plotfiles also contain some minor calculations performed on the final result. The python scrips are called on when using the makefile and finds the necessary datasets automatically.

### `figures/`
This folder contains all the plots used in the report, which are made by the python scripts. It also contains animations of the different experiments we perfomred. 

### `include/`
This folder contains the headerfile `Crank.hpp` for out source file `Crank.hpp`, which is linked in complation. 

### `datafiles/`
This folder is the default save location for all the datafiles created when running the program. Most of the datafiles are bianry files which are read by the pythonscripts using 'pyarma'. As the files are very large, they are not pushed to git, and you will have to run the program yourself if you wan to see them.

### `arkive/`
This folder contains code and other files that was not used in the final version of the program, and can therefore be disregarded.

## Compiling, running and plotting
The program had some high level parallelisation, but might still take a while to run. The datafiles produced by the program are also pretty large, and are thefore not pushed to the repo. To be able to plot you must then first run the program. Below are instuctions on compiling the program, note that you must be in the `project5/` folder for the makefile to run.
### Instuctions:
To compile the program, use the following command:
`$ make`

To run the program, use the following command:
`$ make run `

To generate animations and plots, use the following line:
`$ make plot`

