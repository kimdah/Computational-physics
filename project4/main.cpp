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
  Ising ising(L, T, seed, ordered_spin);
  //Ising ising(20, 10, 12087, 0);
  cout << "The beginning\n";
  ising.print();
  for (int i = 0; i < n_cycles; i++) {
    ising.write_parameters_to_file(ofile);
    for (int j = 0; j < n_cycles/100; j++) {
      ising.run_metropolis_MCMC();
    }
  }
  ofile.close();
  cout << "The end\n";
  ising.print();

  //
  //
  // double totalenergy = energy_of_state(s);
  // double T = 1.0;
  // int n_cycles = 1;
  // int temporary_exp_eps = 0;
  // int temporary_exp_m = 0;
  // int temporary_exp_E = 0;
  // int temporary_exp_M = 0;
  // for(int i=1 ; i<=n_cycles ; i++){
  //   sampling(s, T, totalenergy);
  //   temporary_exp_eps += exp_val_eps_per_cycle;
  //   temporary_exp_m += exp_val_m_per_cycle;
  //   temporary_exp_E +=exp_val_E_per_cycle;
  //   temporary_exp_M += exp_val_M_per_cycle;
  // }
  // double exp_val_eps_all_cycles = temporary_exp_eps/n_cycles;
  // double exp_val_eps_all_cycles_squared = pow(temporary_exp_eps,2)/n_cycles;
  // double exp_val_m_all_cycles = temporary_exp_m/n_cycles;
  // double exp_val_m_all_cycles_squared = pow(temporary_exp_m,2)/n_cycles;
  // double exp_val_E_all_cycles = temporary_exp_E/n_cycles;
  // double exp_val_E_all_cycles_squared = pow(temporary_exp_E,2)/n_cycles;
  // double exp_val_M_all_cycles = temporary_exp_M/n_cycles;
  // double exp_val_M_all_cycles_squared = pow(temporary_exp_M,2)/n_cycles;
  // double headcapacity_all_cycles = (1./N)*(1./pow(T,2))*(exp_val_E_all_cycles_squared
  //    - pow(exp_val_E_all_cycles,2));
  // double susceptibility_all_cycles = (1./N)*(1./pow(T,2))*(exp_val_M_all_cycles_squared
  //    - pow(exp_val_M_all_cycles,2));
  //

  return 0;
}
