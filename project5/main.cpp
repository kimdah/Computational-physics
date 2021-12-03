
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
  
  Crank crankyboii(0.005, 2.5e-4); // randomly chosen!

  //crankyboii.to_file("A");
  //crankyboii.to_file("B");
  crankyboii.to_file("U");

  return 0;
}
