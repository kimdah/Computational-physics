#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h>
#include <vector>
#include <cmath>
#include <iomanip>

#include "./include/Ising.hpp"

using namespace std;
// Performs simulations based on parameter inputs
double simulator(int n_cycles, int lattice_side_length, double T, int seed, int ordered_spin, string filen);
void problem4();
void analytical_2x2(double T);

int main(int argc, char const *argv[]) {
  int T, L, n_cycles, ordered_spin, seed;
  string output_file_name;
   if (argc != 6)
    {
      // Get the name of the executable file
      std::string executable_name = argv[0];

      std::cerr << "Running simulations for Project 4 specified in main.cpp. To run a specific simulation use 5 parameters like so:" << std::endl;
      std::cerr << executable_name << " <Temperature (integer)>"
      << " <lattice side size (integer)>" << " <MCMC cycles (integer)>"
      << " <unordered lattice: use 0, ordered lattice: use -1 or 1>"<< std::endl;
      problem4();
      return 0; // quit program

    } else if (argc == 6) {
      T = atoi(argv[1]);
      L = atoi(argv[2]);
      n_cycles = atoi(argv[3])/100;
      int ordered_spin = atoi(argv[4]); // 0 = unordered, ordered: -1 or 1
      output_file_name = argv[5];
      seed = 2134;
      simulator(n_cycles, L, T, seed, ordered_spin, output_file_name);
    }
  return 0;
}

double simulator(int n_cycles, int lattice_side_length, double T, int seed, int ordered_spin, string output_file_name) {
   // ------ Output-file --------
  string filename = "datafiles/" + output_file_name;
  ofstream ofile;
  ofile.open(filename);
  // Some width and precision parameters we will use to format the output
  int width = 16;
  int prec  = 8;
  ofile << setw(width) << "Sample#";
  ofile << setw(width) << "E";
  ofile << setw(width) << "M";
  ofile << setw(width) << "Expval epsilon";
  ofile << setw(width) << "Expval M";
  ofile << setw(width) << "C_V";
  ofile << setw(width) << "Sucept.";
  ofile << endl;
  // -----------------------------
  // Run the sim
  Ising ising(lattice_side_length, T, seed, ordered_spin);

  int sample_rate = 100;
  // Run MCMC cycles:
  for (int i = 0; i < ((n_cycles/sample_rate)+1); i++) {
    ising.write_parameters_to_file(ofile);
    for (int j = 0; j < sample_rate; j++) {
      ising.run_metropolis_MCMC();
    }
  }
  ofile.close();
  return 42.42;
}

void problem4() {
  // Do all the things we need for Problem 4 here
  int cycles = 20000;
  double temp = 1.0;
  simulator(cycles, 20, temp, 1337, 0, "task4.txt");
  analytical_2x2(temp);
}

void analytical_2x2(double T){  // Maybe in Ising.cpp?
  //double kB = 1.38064852 * pow(10, -23); // Remove?
  double beta = 1. /T;

  double Z = 2*exp(beta*8) + 2*exp(-beta*8) + 12;
  double exp_val_epsilon = (4./Z) * (exp(-beta*8) - exp(beta*8));
  double exp_val_abs_mag = (2./Z) * (exp(beta*8) + 2);
  double heat_capacity = (32./(pow(T,2)*Z))*(exp(-beta*8)-exp(beta*8)- (2./Z)*(exp(-beta*16)+exp(beta*16)-2));
  double susceptibility = (2/(T*Z)) * (exp(8*beta) +1 - ((2./Z)*(exp(16*beta)+ 4*exp(8*beta)+4)));

  // Write to file
  ofstream ofile;
  ofile.open("./datafiles/analytical_2x2_T=" +to_string(T) +".txt");
  int width = 16;
  int prec  = 8;

  ofile << setw(width) << setprecision(prec) << scientific << "T";
  ofile << setw(width) << setprecision(prec) << scientific << "<eps>";
  ofile << setw(width) << setprecision(prec) << scientific << "<m>";
  ofile << setw(width) << setprecision(prec) << scientific << "C_v";
  ofile << setw(width) << setprecision(prec) << scientific << "Chi";
  ofile << endl;

  ofile << setw(width) << setprecision(prec) << scientific << T;
  ofile << setw(width) << setprecision(prec) << scientific << exp_val_epsilon;
  ofile << setw(width) << setprecision(prec) << scientific << exp_val_abs_mag;
  ofile << setw(width) << setprecision(prec) << scientific << heat_capacity;
  ofile << setw(width) << setprecision(prec) << scientific << susceptibility;
  ofile << endl;
  ofile.close();
}
