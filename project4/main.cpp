#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <chrono>
#include "./include/Ising.hpp"

using namespace std;
// Performs simulations based on parameter inputs
double simulator(int n_cycles, int lattice_side_length, double T, int seed, int ordered_spin, string filen);
void problem4();
void problem5_6();
void problem7_8();
void analytical_2x2(double T);
double phase_transitions_parallel(double T_start, double T_end, int steps, int lattice_side_length, int seed, int ordered_spin);
double phase_transitions_serial(double T_start, double T_end, int steps, int lattice_side_length, int seed, int ordered_spin);

template <typename T> string to_string_with_precision(const T a_value, const int n = 1) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

int main(int argc, char const *argv[]) {
  int T, L, n_cycles, ordered_spin, seed;
  string output_file_name;
   if (argc != 6)
    {
      // Get the name of the executable file
      std::string executable_name = argv[0];

      std::cerr << "Running simulations for Project 4 specified in main.cpp. To run a specific simulation use 5 parameters like so:" << std::endl;
      std::cerr << executable_name << " <temperature (float)>"
      << " <lattice side size (integer)>" << " <MCMC cycles (integer)>"
      << " <unordered lattice: use 0, ordered lattice: use -1 or 1>"
      << " <output_file_name> " << std::endl;
      problem4();
      problem5_6();
      //problem7_8();
      return 0; // quit program

    } else if (argc == 6) {
      T = atof(argv[1]);
      L = atoi(argv[2]);
      n_cycles = atoi(argv[3]);
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
  ofile << setw(width) << "<eps>";
  ofile << setw(width) << "<m>";
  ofile << setw(width) << "C_V";
  ofile << setw(width) << "Suscept.";
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
  int cycles = 100000;
  double temp = 1.0;
  simulator(cycles, 2, temp, 1337, 0, "task4.txt"); //unordered
  // Test: simulator(cycles, 40, 2.4, 1337, 0, "task6.txt"); //unordered
  analytical_2x2(temp);
}

void problem5_6() {
  // or use cml arguements
  int cycles = 10000;
  int L = 20;
  double T_1 = 1.0;
  double T_2 = 2.4;
  int seed = 3572;

  // T = 1.0 :
  simulator(cycles, L, T_1, seed, 1, "ncyc_1e5_L_20_T_1.0_ordered.txt"); // for -1 also?
  simulator(cycles, L, T_1, seed, 0, "ncyc_1e5_L_20_T_1.0_unordered.txt");

  // T = 2.4
  simulator(cycles, L, T_2, seed, 1, "ncyc_1e5_L_20_T_2.4_ordered.txt"); // for -1 also?
  simulator(cycles, L, T_2, seed, 0, "ncyc_1e5_L_20_T_2.4_unordered.txt");
}

void problem7_8() {
  // Problem 7: Speedup
  double time_parallel = 0;
  double time_serial = 0;
  auto start = chrono::high_resolution_clock::now();
    //Broad sweeps
  phase_transitions_parallel(1.6, 2.6, 10, 40, 1337, 0);
  auto end = chrono::high_resolution_clock::now();
  chrono::duration<double> p = end - start;
  time_parallel = p.count();

  start = chrono::high_resolution_clock::now();
    //Broad sweeps
  phase_transitions_serial(1.6, 2.6, 10, 40, 1337, 0);
  end = chrono::high_resolution_clock::now();
  chrono::duration<double> s = end - start;
  time_serial = s.count();
  cout << "Parallel time: " << time_parallel << " Serial time: " << time_serial << "\n";
  cout << "Speedup = " << time_serial/time_parallel << "\n";
  // Problem 8: Critical T
  //Broad sweeps of T=2.1 to T=2.4
  //L=40
  phase_transitions_parallel(2.1, 2.4, 10, 40, 1337, 0);
  //L=60
  phase_transitions_parallel(2.1, 2.4, 10, 60, 1337, 0);
  //L=80
  phase_transitions_parallel(2.1, 2.4, 10, 80, 1337, 0);
  //L=100
  phase_transitions_parallel(2.1, 2.4, 10, 100, 1337, 0);
}


void analytical_2x2(double T){
  //double kB = 1.38064852 * pow(10, -23); // Remove?
  double beta = 1./T;

  double Z = 2*exp(beta*8) + 2*exp(-beta*8) + 12;
  double exp_val_epsilon = (4./Z) * (exp(-beta*8) - exp(beta*8));
  double exp_val_abs_mag = (2./Z) * (exp(beta*8) + 2);
  // test
  double eps2 = (8./Z) * (exp(-beta*8) + exp(beta*8));
  double m2 = (1./Z) * (2*exp(beta*8) + 2);
  //
  double heat_capacity = (32./(pow(T,2)*Z))*(exp(-beta*8)+exp(beta*8) - ((2./Z)*(exp(-beta*16)+exp(beta*16)-2)));
  double susceptibility = (8./(T*Z)) * (exp(8*beta) + 1 - ((2./Z)*(exp(16*beta)+ 4*exp(8*beta)+4)));

  // Write to file
  ofstream ofile;
  //ofile.open("./datafiles/analytical_2x2_T=" +to_string(T) +".txt");

  // To have 2 decimals in output-filename
  std::ostringstream temp;
  temp << std::fixed << std::setprecision(1) << T;

  ofile.open("./datafiles/analytical_2x2_T=" + temp.str() +".txt");// to_string(T)
  int width = 18;
  int prec  = 8;

  ofile << setw(2)<< "T" << setw(width) << "<eps>"
  << setw(width) << "<m>" << setw(width) << "C_v"<< setw(width) << "chi" << endl;

  ofile << setprecision(2) << scientific << T;
  //ofile << setw(width) << setprecision(prec) << scientific << Z;
  ofile << setw(width) << setprecision(prec) << scientific << exp_val_epsilon;
  //ofile << setw(width) << setprecision(prec) << scientific << eps2;
  ofile << setw(width) << setprecision(prec) << scientific << exp_val_abs_mag;
  //ofile << setw(width) << setprecision(prec) << scientific << m2;
  ofile << setw(width) << setprecision(prec) << scientific << heat_capacity;
  ofile << setw(width) << setprecision(prec) << scientific << susceptibility;
  ofile << endl;
  ofile.close();
}

double phase_transitions_parallel(double T_start, double T_end, int steps, int lattice_side_length, int seed, int ordered_spin){
  double h = (T_end - T_start) / steps;
  int burn_in = 10000;
  int samples = 100;
  double time = 0;
  //cout << "starting parallel region \n";
  string filename = "datafiles/phase_transitions_parallel_T("+to_string_with_precision(T_start)+"-"+to_string_with_precision(T_end)+")_L("+to_string(lattice_side_length)+")_steps("+to_string(steps)+").txt";
  ofstream ofile;
  ofile.open(filename);
  int width = 16;
  ofile << setw(width) << "T";
  ofile << setw(width) << "<eps>";
  ofile << setw(width) << "<m>";
  ofile << setw(width) << "C_V";
  ofile << setw(width) << "Sucept.";
  ofile << endl;

  #pragma omp parallel for
  for (int i = 0; i < steps+1; i++){
    //cout << "Process doing stuff " << i << "\n";
    double T = T_start + i * h;
    Ising ising(lattice_side_length, T, seed, ordered_spin);
    for (int j = 0; j < burn_in; j++) {
      ising.run_metropolis_MCMC();
    }
    #pragma omp critical
    ising.sample_average_over_sampled_values(ofile, samples);
  }
  ofile.close();
  return time;
}

double phase_transitions_serial(double T_start, double T_end, int steps, int lattice_side_length, int seed, int ordered_spin) {
  double h = (T_end - T_start) / steps;
  int burn_in = 10000;
  int samples = 100;
  double time = 0;
  //cout << "starting parallel region \n";
  string filename = "datafiles/ignore.txt";
  ofstream ofile;
  ofile.open(filename);
  int width = 16;
  ofile << setw(width) << "T";
  ofile << setw(width) << "<eps>";
  ofile << setw(width) << "<m>";
  ofile << setw(width) << "C_V";
  ofile << setw(width) << "Sucept.";
  ofile << endl;
  for (int i = 0; i < steps+1; i++){
    //cout << "Process doing stuff " << i << "\n";
    double T = T_start + i * h;
    Ising ising(lattice_side_length, T, seed, ordered_spin);
    for (int j = 0; j < burn_in; j++) {
      ising.run_metropolis_MCMC();
    }
    ising.sample_average_over_sampled_values(ofile, samples);
  }
  ofile.close();
  return time;
}
