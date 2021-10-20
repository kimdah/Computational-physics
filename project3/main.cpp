#include <string>
#include <iostream>
#include <armadillo>
#include <string>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <sstream>



//Code files included:
#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace std;
using namespace arma;

int main(int argc, char const *argv[]) {

  double d = pow(10,4);

  PenningTrap penning_trap(9.65*10, 9.65*pow(10,8), pow(10,4));

  // double q_in, double m_in, arma::vec pos_in, arma::vec vel_in
  //Particle new_particle(1, 40.08, vec(3).randn()*d*0.1, vec(3).randn()*d*0.1); // Ca ATOM!

  //penning_trap.add_particle(new_particle);


  return 0;
}
