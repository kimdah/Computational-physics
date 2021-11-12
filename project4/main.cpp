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
  
  vector< vector<int> > s_candid;
  for (int i = 0; i < L; i++){
    vector<int> s_row(L, 1);
    s_candid.push_back(s_row);
  }

  // random from pr3: vec(3).randn()*0.1*penning_trap.d_
  int randRow = rand() % L;
  int randCol = rand() % L;
  s_candid[randRow][randCol] *= -1;

  std::cout << s_candid[randRow][randCol]<< endl;
  





  // for (int i = 0; i < L; i++){
  //   for (int j = 0; j < L; j++){
  //     //s_config()
  //   }
  // }



  return 0;
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
