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

  vector< vector<int> > s ; // sample
  for (int i = 0; i < L; i++){
    vector<int> s_row(L, 1);
    s.push_back(s_row);
  }

  double totalenergy = energy_of_state(s);
  double T = 1.0;
  int n_cycles = 1;
  for(int i=1 ; i<=n_cycles ; i++){
    sampling(s, T, totalenergy);
  }



  return 0;
}

std::vector<std::vector<int> > sampling(std::vector<std::vector<int> > s_current, double T, double totalenergy){
  // running one MC cycle for sampling
  std::vector<double> boltzmann_factors = boltzmann_factor(T);

  for (int c = 0; c < N; c++){ // one MC cycle
    // flip random spin
    int randRow = rand() % L;
    int randCol = rand() % L;
    s_current[randRow][randCol] *= -1;

    // find p(s')/p(s_i)=e^(-beta(E(sdash)-(E(s_i))=e^(-beta*deltaE)

    // examining surrounding spins to figure out index in boltzmann_factor vector
    // for computing the probability ratio

    int sumofsurroundingspins = s_current[(randRow - 1 + L) % L][randCol] // Neighbour below
                              + s_current[randRow][(randCol - 1 + L) % L] // Neighbour to the left
                              + s_current[(randRow + 1) % L][randCol] // Neighbour above
                              + s_current[randRow][(randCol + 1) % L]; // Neighbour to the right


    
    // finding the index to use in Boltzmann
    int index;
    int deltaE = sumofsurroundingspins*2;
    //int m = sumofsurroundingspins;

    // boltzmann factor depends on flipping a +1 to -1, so the value will have
    // reverse index when a negative spin is flipped to positive.
    if(s_current[randRow][randCol] == 1){
      index = 5 - sumofsurroundingspins/2 + 2; // reversing index
    }else{
      index = sumofsurroundingspins/2 + 2;
    }

    // Acceptance ratio
    double probability_ratio = boltzmann_factors[index];
    double r = rand()/RAND_MAX;

    if(r > probability_ratio){ //Rejected spin-flip
      s_current[randRow][randCol] *= -1;
    }
    else{
      // Accept spin configuration candidate
      double totalenergy = totalenergy + deltaE;
      double epsilon = totalenergy/N;
    }

    /*
    double randRow_index = (randRow + L)%L;
    double randCol_index = (randCol + L)%L;
    double sumofsurroundingspins = s_current[randRow_index-1][randCol] + s_current[randRow][randCol_index-1] + s_current[randRow_index+1][randCol] + s_currenct[randRow]
    [randCol_index+1];
    */
    /* // Old code for summing neighbour spins
    int sumofsurroundingspins; // NOT energy, simply sum to get a sense of how many +1 neighbour spins there are
    if(randRow != 0 && randRow != L-1 && randCol != 0 && randCol != L-1){
      sumofsurroundingspins = s_current[randRow-1][randCol] + s_current[randRow][randCol-1]
       + s_current[randRow+1][randCol] + s_current[randRow][randCol+1];
    }
    else if(randRow==0 && randCol==0){
      sumofsurroundingspins = s_current[L-1][randCol] + s_current[randRow][L-1]
       + s_current[randRow+1][randCol] + s_current[randRow][randCol+1];
    }
    else if(randRow==L-1 && randCol==L-1){
      sumofsurroundingspins = s_current[randRow-1][randCol] + s_current[randRow][randCol-1]
       + s_current[0][randCol] + s_current[randRow][0];
    }
    else if(randRow==0 && randCol==L-1){
      sumofsurroundingspins = s_current[L-1][randCol] + s_current[randRow][randCol-1]
       + s_current[randRow+1][randCol] + s_current[randRow][0];
    }
    else if(randRow==L-1 && randCol==0){
      sumofsurroundingspins = s_current[randRow-1][randCol] + s_current[randRow][L-1]
       + s_current[0][randCol] + s_current[randRow][randCol+1];
    }
    else if(randRow==0){
      sumofsurroundingspins = s_current[L-1][randCol] + s_current[randRow][randCol-1]
       + s_current[randRow+1][randCol] + s_current[randRow][randCol+1];
    }
    else if(randCol==0){
      sumofsurroundingspins = s_current[randRow-1][randCol] + s_current[randRow][L-1]
       + s_current[randRow+1][randCol] + s_current[randRow][randCol+1];
    }
    else if(randRow==L-1){
      sumofsurroundingspins = s_current[randRow-1][randCol] + s_current[randRow][randCol-1]
       + s_current[0][randCol] + s_current[randRow][randCol+1];
    }
    else if(randCol==L-1){
      sumofsurroundingspins = s_current[randRow-1][randCol] + s_current[randRow][randCol-1]
       + s_current[randRow+1][randCol] + s_current[randRow][0];
    }
 */

  }
  return s_current;
}

double energy_of_state(std::vector<std::vector<int> > s){
  // finding the energy of a particular spin configuration s
  double energy;
  for(int i=1 ; i<L+1 ; i++){ //the first row will be the Lth row
    for(int j=1 ; j<L+1 ; j++){ //the first column will be the Lth column
      int i_index = (i + L)%L;
      int j_index = (j + L)%L;
      energy += s[i_index][j_index]*s[i_index-1][j_index] + s[i_index][j_index]*s[i_index][j_index-1];
    }
  }
  return energy;
}



std::vector<double> boltzmann_factor(double T){
  double beta = 1. / T;
  vector<double> boltzmann_values;
  boltzmann_values.push_back(exp(-beta*(-8))); // 0 +1 spins
  boltzmann_values.push_back(exp(-beta*(-4))); // 1 +1 spins
  boltzmann_values.push_back(exp(-beta*(0))); // = 1? // 2 +1 spins
  boltzmann_values.push_back(exp(-beta*(4))); // 3 +1 spins
  boltzmann_values.push_back(exp(-beta*(8))); // 4 +1 spins
  return boltzmann_values;
}
