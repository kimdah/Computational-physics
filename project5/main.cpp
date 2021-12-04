
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
  crank.to_file("U");
  cx_cube nicolson = crank.run_simulation(321);
  crank.output_probabilities(nicolson, "datafiles/probability_sum_test.txt");
}
