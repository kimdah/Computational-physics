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





std::vector<double> boltzmann_factor(double T){
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
  vector<double> boltzmann_factors = boltzmann_factor(temp);
  double Z = 2*boltzmann_factors[-]

}



analytical_energy_per_spin_2x2(){

}
