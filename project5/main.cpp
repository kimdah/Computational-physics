
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>
#include "./include/Crank.hpp"
#include <armadillo>

#include <complex>
#include <cmath>

using namespace std::complex_literals; // to use imaginary number i | DEMANDS c++14!
using namespace std;
using namespace arma;

void problem7();
int main(int argc, char const *argv[]) {
  // Crank(double h, double deltat)

  problem7();

  return 0;
}

void problem7() {
  Crank crank(0.005, 2.5e-5); // NOT randomly chosen

  //crank.to_file("A");
  //crank.to_file("B");
  
  cx_cube prob7_1 = crank.run_simulation(321); // TODO: 321 gives T = 0.008. Change input to actual time and modify code in Crank
  crank.to_file("U");
  crank.output_probabilities(prob7_1, "datafiles/probability_sum_test.txt");
  prob7_1.save("datafiles/prob7_1");
}
