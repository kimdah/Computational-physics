#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <cmath>
#include <iomanip>

#include "./include/Ising.hpp"

using namespace std;


int main(int argc, char const *argv[]) {
   if (argc != 6)
    {
      // Get the name of the executable file
      std::string executable_name = argv[0];

      std::cerr << "Error: Wrong number of input arguments. 5 expected." << std::endl;
      std::cerr << "Usage: " << executable_name << " <Temperature (integer)>"
      << " <lattice side size (integer)>" << " <MCMC cycles (integer)>"
      << " <unordered lattice: use 0, ordered lattice: use -1 or 1>"<< std::endl;
      return 1;
    }

  // Read command line arguments
  const int T = atoi(argv[1]);
  const int L = atoi(argv[2]);
  const int n_cycles = atoi(argv[3])/100;
  const int ordered_spin = atoi(argv[4]); // 0 = unordered, ordered: -1 or 1
  const string output_file_name = argv[5];
  const int seed = 2134;

  // ------ Output-file --------
  string filename = "datafiles/" + output_file_name;
  ofstream ofile;
  ofile.open(filename);
  // Some width and precision parameters we will use to format the output
  int width = 16;
  int prec  = 8;
  ofile << setw(width) << setprecision(prec) << scientific << "Sample#";
  ofile << setw(width) << setprecision(prec) << scientific << "E";
  ofile << setw(width) << setprecision(prec) << scientific << "M";
  ofile << endl;
  // -----------------------------

  Ising ising(L, T, seed, ordered_spin);
  //Ising ising(20, 10, 12087, 0);
  cout << "The beginning\n";
  ising.print();
  // Run MCMC cycles:
  for (int i = 0; i < n_cycles; i++) {
    ising.write_parameters_to_file(ofile);
    for (int j = 0; j < n_cycles/100; j++) {
      ising.run_metropolis_MCMC();
      // get epsilon for each cycle here
      //cout << "<eps> for " << i << " cycles" << ising.expval_eps(i) << endl;

    }
  }
  cout << "<eps>" << ising.expval_eps(n_cycles) << endl;
  ofile.close();
  cout << "The end\n";
  ising.print();

  return 0;
}

//-------- TRIED TO GET THE ANALYTICAL VALUES AND ESTIMATES FOR 2X2 TO FILE - NOT WORKING

// void test_with_2x2(n_cycles){
//   ofstream tfile;
//   tfile.open("values_table_2x2.txt");
//   int width = 16;
//   int prec  = 8;
//   tfile << setw(width) << setprecision(prec) << scientific << "Sample#";
//   tfile << setw(width) << setprecision(prec) << scientific << "E";
//   tfile << setw(width) << setprecision(prec) << scientific << "M";
//   tfile << endl;
//   analytical_2x2(tfile);
//
//   Ising ising(2, 1.0, 2134, 0);
//
//   for (int i = 0; i < n_cycles; i++) {
//     ising.analytical_2x2(tfile);
//     // for (int j = 0; j < n_cycles/100; j++) {
//     //   ising.run_metropolis_MCMC();
//     // }
//   }
//   tfile.close();
// }
//
// void analytical_2x2(){  // Maybe in Ising.cpp?
//   vector<double> analytical_values;
//   double kB = 1.38064852 * pow(10, -23); // Remove?
//   double beta = 1. / (kB*T_);
//
//   double Z = 2*exp(beta*8) + 2*exp(-beta*8) + 12;
//   double exp_val_epsilon = (4./Z) * (exp(-beta*8) - exp(beta*8));
//   double exp_val_abs_mag = (2./Z) * (exp(beta*8) + 2);
//   double heat_capacity = (32./(kB*pow(T,2)*Z))*(exp(-beta*8)-exp(beta*8)- (2./Z)*(exp(-beta*16)+exp(beta*16)-2));
//   double susceptibility = (2/(kB*T*Z)) * (exp(8*beta) +1 - ((2./Z)*(exp(16*beta)+ 4*exp(8*beta)+4)));
// }
