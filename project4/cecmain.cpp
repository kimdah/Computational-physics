#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <cmath>

using namespace std;

std::vector<double> boltzmann_factor(double T);

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
  int L = 2;
  int N = L*L;

  vector< vector<int> > s_current ; // sample
  for (int i = 0; i < L; i++){
    vector<int> s_row(L, 1);
    s_current.push_back(s_row);
  }




  //std::cout << s_candidate[randRow][randCol] << endl;




  return 0;
}

std::vector<double> sampling(std::vector<double> s){
  for (int c = 0; c < N; c++){ // one MC cycle
    // flip random spin
    int randRow = rand() % L;
    int randCol = rand() % L;
    s_current[randRow][randCol] *= -1;

    // find p(s')/p(s_i)
  }
}

energy_of_state(std::vector<double> s){

}



std::vector<double> boltzmann_factor(double T){
  // possible boltzmann factors for the 5 possible energy differences when flipping spins
  double kB = 1.38064852 * pow(10, -23); //
  double beta = 1. / (kB*T);
  vector<double> boltzmann_values;
  boltzmann_values.push_back(exp(-beta*(-8))); // 0 +1 spins
  boltzmann_values.push_back(exp(-beta*(-4))); // 1 +1 spins
  boltzmann_values.push_back(exp(-beta*(0))); // = 1? // 2 +1 spins
  boltzmann_values.push_back(exp(-beta*(4))); // 3 +1 spins
  boltzmann_values.push_back(exp(-beta*(8))); // 4 +1 spins
  return boltzmann_values;

}

analytical_2x2(double temp){
  double kB = 1.38064852 * pow(10, -23); //
  double beta = 1. / (kB*T);


  double Z = 2*exp(beta*8) + 2*exp(-beta*8) + 12;
  double exp_val_epsilon = (4./Z) * (exp(-beta*8) - exp(beta*8));
  double exp_val_abs_mag = (2./Z) * (exp(beta*8) + 2)
  double heat_capacity = (32./(kB*pow(temp,2)*Z))*(exp(-beta*8)-exp(beta*8)- (2./Z)*(exp(-beta*16)+exp(beta*16)-2));


}

void expval_energy_per_spin(){

}

void Ising::alternative_sampling(std::vector<std::vector<int> > s_current, double T){
  // cecilie
  std::vector<double> boltzmann_factors = boltzmann_factor(T);

  double energy_before = calc_tot_energy_of_state(s_current);

  for (int i = 0; i < N; i++){ // one MC cycle
    // flip random spin
    int randRow = rand() % L;
    int randCol = rand() % L;

    s_current[randRow][randCol] *= -1; // flip spin

    double energy_after = calc_tot_energy_of_state(s_current);
    double delta_energy = energy_after - energy_before;
    // ....
  }
}
