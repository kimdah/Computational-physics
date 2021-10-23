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

void problem_10(double f, vec frequencies, bool interactions);
double simulator(int iterations, int duration, int particles, std::string outputs, bool interactions, bool euler_cromer, bool pertrubation, bool randomseed, double f, double w_v, bool out);

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 1)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

int main(int argc, char const *argv[]) {

    // ------------- PROBLEM 9 -----------
    //Problem 9 point 1
    simulator(10000, 100, 1, "tz", true, false, false, false, 0.0, 0.0, true);
    // Problem 9 point 2
    // interactions on
    simulator(10000, 100, 2, "xy", true, false, false, false, 0.0, 0.0, true);
    // interactions off
    simulator(10000, 100, 2, "xy", false, false, false, false, 0.0, 0.0, true);
    // Point 3: Phase space plots
    simulator(10000, 100, 2, "xv", true, false, false, false, 0.0, 0.0, true);
    simulator(10000, 100, 2, "xv", false, false, false, false, 0.0, 0.0, true);
    simulator(10000, 100, 2, "yv", true, false, false, false, 0.0, 0.0, true);
    simulator(10000, 100, 2, "yv", false, false, false, false, 0.0, 0.0, true);
    simulator(10000, 100, 2, "zv", true, false, false, false, 0.0, 0.0, true);
    simulator(10000, 100, 2, "zv", false, false, false, false, 0.0, 0.0, true);
    // Point 4: 3D plot
    simulator(10000, 100, 2, "xyz", true, false, false, false, 0.0, 0.0, true);
    simulator(10000, 100, 2, "xyz", false, false, false, false, 0.0, 0.0, true);
    // Point 5: step sizes
    for (int i = 1; i < 6; i++) {
        simulator(pow(10,i), 100, 1, "txyz", true, false, false, false, 0.0, 0.0, true); // RK4
        simulator(pow(10,i), 100, 1, "txyz", true, true, false, false, 0.0, 0.0, true); // Euler Cromer
    }

    // ------------- PROBLEM 10 -----------
    // For each of the amplitudes f=0.1,0.4,0.7, produce a graph that shows the fraction of
    // particles that are still trapped after 500μs as a function of the applied angular frequency ω_V


    vector<double> freqs1;
    for (double i = 0.2; i<2.51; i+=0.05) {
        freqs1.push_back(i); //*pow(10,6)
    }
    
    
    problem_10(0.4, freqs1, false);
    problem_10(0.7, freqs1, false);
    problem_10(0.1, freqs1, false);
    vector<double> freqs2;
    for (double i = 0.2; i<0.81; i+=0.001) {
        freqs2.push_back(i); //*pow(10,6)
    }
    problem_10(0.1, freqs2, false);
    problem_10(0.4, freqs2, false);
    problem_10(0.7, freqs2, false);
    return 0;
}

void problem_10(double f, vec frequencies, bool interactions) {
    int width = 16;
    int prec  = 8;
    std::string filename = "Results/problem10_f_"+to_string_with_precision(f)+".txt";
    std::ofstream ofile;
    ofile.open(filename);
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << "w_v";
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << "fraction";
    ofile << std::endl;
    for (double x : frequencies){
        std::cout << x << "----------------------" << std::endl;
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x;
        double fraction = simulator(5000, 500, 100, "w", interactions, false, true, true, f, x, false);
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << fraction;
        ofile << std::endl;
    }
    std::cout << std::endl;
    ofile.close();

}
// Simulator with parameters in order the number of iterations, duration, particles, outputs included, whether particle interactions are on, EC on/off,
// V-pertrubation for Problem 10, random seed on/off, V-pertrubation amplitude, V-pertrubation angular frequency
// Returns the fraction of particles still trapped at the end of the simulation
double simulator(int iterations, int duration, int particles, std::string outputs, bool interactions, bool euler_cromer, bool pertrubation, bool randomseed, double f, double w_v, bool out) {

    int n = iterations; //itrs
    //std::cout << "iterations " << n;
    double t = duration;
    double h = double(t) / n;
    double d = pow(10,4);
    // Create a filename specifying in order the number of iterations, duration, particles, outputs included, whether particle interactions are on, EC on/off,
    // V-pertrubation for Problem 10, random seed on/off, V-pertrubation amplitude, V-pertrubation angular frequency
    std::string filename;
    if (euler_cromer && out) {
            filename = "Results/EC_i_"+std::to_string(iterations)+"_d_"+std::to_string(duration)+"_p_"+std::to_string(particles)+"_pi_"+std::to_string(interactions)+
            "_outputs_"+outputs+"_pert_"+std::to_string(pertrubation)+"_rs_"+std::to_string(randomseed)+"_f_"+to_string_with_precision(f)+"_w_v_"+to_string_with_precision(w_v)+".txt";
        } else if (out){
            filename = "Results/RK4_i_"+std::to_string(iterations)+"_d_"+std::to_string(duration)+"_p_"+std::to_string(particles)+"_pi_"+std::to_string(interactions)+
            "_outputs_"+outputs+"_pert_"+std::to_string(pertrubation)+"_rs_"+std::to_string(randomseed)+"_f_"+to_string_with_precision(f)+"_w_v_"+to_string_with_precision(w_v)+".txt";
        }

    std::ofstream ofile;
    ofile.open(filename);

    // Some width and precision parameters we will use to format the output
    int width = 16;
    int prec  = 8;


    // Run simulation for a single particle with default values for B_0, V_0 and d
    // double B0_in, double V0_in, double d_in
    PenningTrap penning_trap(9.65*10, 9.65*pow(10,8), pow(10,4));

    // For problem 10:
    if (pertrubation) {
        penning_trap.V0_= 2.412131395*pow(10,5);
        penning_trap.E_= 2.412131395*pow(10,5);
        penning_trap.d_ = 500;
    }

    // Set random seed so that we have comparable results
    
    arma_rng::set_seed_random();
    // Turn on or off particle/coloumb interactions
    penning_trap.particle_interactions_ = interactions;
    // double q_in, double m_in, arma::vec pos_in, arma::vec vel_in
    for (int j = 0; j < particles; j++) {
        if(!randomseed) {arma_rng::set_seed(j*3);}
        Particle new_particle(1, 40.078, vec(3).randn()*0.1*penning_trap.d_, vec(3).randn()*0.1*penning_trap.d_); // Ca ATOM!
        penning_trap.add_particle(new_particle);
    }


    if (outputs.find('t') != std::string::npos) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << "t";}

    for (int j = 1; j < penning_trap.particles_.size()+1; j++) {
        if (outputs.find('x') != std::string::npos && out) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << "x_" + std::to_string(j);}
        if (outputs.find('y') != std::string::npos && out) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << "y_" + std::to_string(j);}
        if (outputs.find('z') != std::string::npos && out) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << "z_" + std::to_string(j);}
        if (outputs.find('v') != std::string::npos && outputs.find('x') != std::string::npos && out) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << "v_x" + std::to_string(j);}
        if (outputs.find('v') != std::string::npos && outputs.find('y') != std::string::npos && out) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << "v_y" + std::to_string(j);}
        if (outputs.find('v') != std::string::npos && outputs.find('z') != std::string::npos && out) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << "v_z" + std::to_string(j);}
    }
    ofile << std::endl;
    if (outputs.find('t') != std::string::npos && out) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << 0;}

    for (int j = 0; j < penning_trap.particles_.size(); j++) {
        if (outputs.find('x') != std::string::npos && out) {ofile << std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].pos_[0];} //endret particles_0 =particles_j
        if (outputs.find('y') != std::string::npos && out) {ofile << std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].pos_[1];}
        if (outputs.find('z') != std::string::npos && out) {ofile << std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].pos_[2];}
        if (outputs.find('v') != std::string::npos && outputs.find('x') != std::string::npos && out) {ofile << std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].vel_[0];}
        if (outputs.find('v') != std::string::npos && outputs.find('y') != std::string::npos && out) {ofile << std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].vel_[1];}
        if (outputs.find('v') != std::string::npos && outputs.find('z') != std::string::npos && out) {ofile << std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].vel_[2];}

    }
    ofile<< std::endl;

    for (int i = 1; i < n+1; i++) {
        if (euler_cromer) {
            penning_trap.evolve_Euler_Cromer(h);
        } else {
            penning_trap.evolve_RK4(h);
        }

        if (pertrubation) {
            //std::cout << "pertrubating " << penning_trap.V0_;
            penning_trap.pertrubation = true;
            penning_trap.E_ = penning_trap.V0_ * (1 + f * cos (w_v * h*i)); //*pow(10,6)
            //std::cout << penning_trap.V0_;
        }

        if (outputs.find('t') != std::string::npos) {ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << h*i;}
        for (int j = 0; j < penning_trap.particles_.size(); j++) {
            if (outputs.find('x') != std::string::npos && out) {ofile << std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].pos_[0];}
            if (outputs.find('y') != std::string::npos && out) {ofile << std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].pos_[1];}
            if (outputs.find('z') != std::string::npos && out) {ofile << std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].pos_[2];}
            if (outputs.find('v') != std::string::npos && outputs.find('x') != std::string::npos && out) {ofile << std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].vel_[0];}
            if (outputs.find('v') != std::string::npos && outputs.find('y') != std::string::npos && out) {ofile << std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].vel_[0];}
            if (outputs.find('v') != std::string::npos && outputs.find('z') != std::string::npos && out) {ofile << std::setw(width) << std::setprecision(prec) << std::scientific << penning_trap.particles_[j].vel_[0];}
        }
        ofile << std::endl;
      }

      ofile.close();
      //std::cout << "extreme " << penning_trap.extreme_ << std::endl;
      return penning_trap.particles_inside()/penning_trap.particles_.size();
}
