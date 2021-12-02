
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
  // Crank(int M, double h, double deltat, double r, double v0)
  Crank crankyboii(5, 0.1, 0.1, 2, 0); // randomly chosen!

  crankyboii.print();
}
