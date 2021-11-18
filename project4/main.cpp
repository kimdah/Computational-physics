#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <cmath>

#include "./include/Ising.hpp"

using namespace std;



int main(int argc, char const *argv[]) {

  // ------ FROM EXAMPLE CODE-------
  // Check number of command line arguments
  assert(argc == 5);

  // // Read command line arguments
  const int L = atoi(argv[1]);
  const int n_cycles = atoi(argv[2]);
  const int ordered_spin = atoi(argv[3]); // 0 = unordered, ordereded: -1 or 1
  const string output_file_name = argv[4];
  // Prepare for file output
  // const int print_prec = 10;
  // ofstream ofile;
  // ofile.open(output_file_name.c_str(), ofstream::trunc);


  // Monte Carlo temporary
  // L = 2;
  int N = L*L;
  int T = 1.0;
  int seed = 1; // ?
  Ising ising(L, T, seed, ordered_spin);




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
