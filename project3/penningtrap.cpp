#include "PenningTrap.hpp";



// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in);

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in);

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r);

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r);

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j);

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i);

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i);

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i);

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt);

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt);
