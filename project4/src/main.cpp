#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <cmath>

using namespace std;

std::vector<double> boltzmann_factor(double T);
std::vector<std::vector<int> > sampling(std::vector<std::vector<int> > s_current, double T, double totalenergy);
double energy_of_state(std::vector<std::vector<int> > s);

int L;
int N;

int main(int argc, char const *argv[]) {

  // ------ FROM EXAMPLE CODE-------
  // Check number of command line arguments
  // assert(argc == 4);

  // // Read command line arguments
  // const int n_A = atoi(argv[1]);
  // const int n_cycles = atoi(argv[2]);
  // const string output_file_name = argv[3];
  // Prepare for file output
  const int print_prec = 10;
  // ofstream ofile;
  // ofile.open(output_file_name.c_str(), ofstream::trunc);


  // Monte Carlo temporary
  L = 2;
  N = L*L;



  double totalenergy = energy_of_state(s);
  double T = 1.0;
  int n_cycles = 1;
  int temporary_exp_eps = 0;
  for(int i=1 ; i<=n_cycles ; i++){
    sampling(s, T, totalenergy);
    temporary_exp_eps += exp_val_eps_per_cycle;
  }
  int exp_val_eps_all_cycles = temporary_exp_eps/n_cycles;



  return 0;
}
