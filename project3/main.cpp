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

    int n = 100; //itrs
    double t = 100.0;
    double h = t / n;
    double d = pow(10,4);
    std::string filename = "Results/9_1_single_particle.txt";
    std::ofstream ofile;
    ofile.open(filename);

    // Some width and precision parameters we will use to format the output
    int width = 16;
    int prec  = 8;

    // Run simulation for a single particle with default values for B_0, V_0 and d
    // double B0_in, double V0_in, double d_in
    PenningTrap penning_trap(9.65*10, 9.65*pow(10,8), pow(10,4));

    // double q_in, double m_in, arma::vec pos_in, arma::vec vel_in
    Particle new_particle(1, 40.08, vec(3).randn()*d*0.1, vec(3).randn()*d*0.1); // Ca ATOM!

    penning_trap.add_particle(new_particle);

    for (int i = 0; i < n+1; i++) {
        penning_trap.evolve_RK4(h);
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << h*i
              << std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[0].pos_[2]
              << std::endl;
    }
    return 0;
}
