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

void simulator(int iterations, int duration, int particles, std::string axis, bool interactions);

int main(int argc, char const *argv[]) {

    // ------------- PROBLEM 9 -----------
    //Problem 9 point 1
    simulator(100, 100, 1, "z", true);
    // Problem 9 point 2
    // interactions on
    simulator(100, 100, 2, "xy", true);
    // interactions off
    simulator(100, 100, 2, "xy", false);



    return 0;
}

void simulator(int iterations, int duration, int particles, std::string axis, bool interactions) {
    const auto test = axis;
    int n = iterations; //itrs
    double t = duration;
    double h = t / n;
    double d = pow(10,4);
    // Create a filename specifying in order the number of iterations, duration, particles, axes included, whether particle interactions are on
    std::string filename = "Results/i_"+std::to_string(iterations)+"_d_"+std::to_string(duration)+"_p_"+std::to_string(particles)+"_pi_"+std::to_string(interactions)+"_axis_"+axis+".txt";
    std::ofstream ofile;
    ofile.open(filename);

    // Some width and precision parameters we will use to format the output
    int width = 16;
    int prec  = 8;


    // Run simulation for a single particle with default values for B_0, V_0 and d
    // double B0_in, double V0_in, double d_in
    PenningTrap penning_trap(9.65*10, 9.65*pow(10,8), pow(10,4));
    penning_trap.particle_interactions_ = interactions;
    // double q_in, double m_in, arma::vec pos_in, arma::vec vel_in
    for (int j = 0; j<particles; j++) {
        Particle new_particle(1, 40.08, vec(3).randn()*d*0.1, vec(3).randn()*d*0.1); // Ca ATOM!
        penning_trap.add_particle(new_particle);
    }

    ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << "t";

    for (int j = 1; j < penning_trap.particles_.size()+1; j++) {
        if (axis.find('x') != std::string::npos) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << "x_" + std::to_string(j);}
        if (axis.find('y') != std::string::npos) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << "y_" + std::to_string(j);}
        if (axis.find('z') != std::string::npos) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << "z_" + std::to_string(j);}
    }
    ofile << std::endl;
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << 0;
    for (int j = 0; j < penning_trap.particles_.size(); j++) {
        if (axis.find('x') != std::string::npos) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[0].pos_[0];}
        if (axis.find('y') != std::string::npos) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[0].pos_[1];}
        if (axis.find('z') != std::string::npos) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[0].pos_[2];}
    }

    ofile<< std::endl;
    for (int i = 1; i < n+1; i++) {
        penning_trap.evolve_RK4(h);
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << h*i;
        for (int j = 0; j < penning_trap.particles_.size(); j++) {
            if (axis.find('x') != std::string::npos) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].pos_[0];}
            if (axis.find('y') != std::string::npos) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].pos_[1];}
            if (axis.find('z') != std::string::npos) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].pos_[2];}
        }
        ofile << std::endl;
      }
    ofile.close();
}
