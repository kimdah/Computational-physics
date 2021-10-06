#include "PenningTrap.hpp"
#include "Particle.hpp"


//q = particle_[0].q_;


// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in){
  B0_ = B0_in; // definer disse
  V0_ = V0_in;
  d_ = d_in;
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in){
  particles_.push_back(p_in);
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r){
  arma::vec E_field = arma::vec(3).fill(0);
  E_field(0) = (r(0)*V0)/((d_)^2);
  E_field(1) = (r(1)*V0)/((d_)^2);
  E_field(2) = (-2*r(2)*V0)/((d_)^2);
  return arma::vec(1);

}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r){
  return arma::vec(1);

}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j){



  return arma::vec(1);

}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i){
  return arma::vec(1);

}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i){
  return arma::vec(1);

}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i){
  return arma::vec(1);
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt){

}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt){

}
