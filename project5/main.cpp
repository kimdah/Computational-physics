
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


int main(int argc, char const *argv[]) {
  // Crank(double h, double deltat)

  Crank crank(0.005, 2.5e-4); // randomly chosen!

  //crank.to_file("A");
  //crank.to_file("B");
  crank.to_file("U");

  return 0;
}
