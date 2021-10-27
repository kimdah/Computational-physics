#include "Particle.hpp"


// Constructor
Particle::Particle(double q_in, double m_in, arma::vec pos_in, arma::vec vel_in)
{
  q_ = q_in;          // charge of particle, [e] elementary charge
  m_ = m_in;          // mass of particle in atomic mass unit u
  pos_.swap(pos_in);  // position vector r
  vel_.swap(vel_in);  // velocity vector v
  outofbounds_ = false; // Is the particle out of bounds?

}

/*
If you cant it to be private

// Method that returns the charge
int Particle::charge(){
  return q_;
}

// Method that returns the mass_
double Particle::mass(){
  return m_;
}

// Method that returns the position vector of the particle
arma::vec Particle::position(){
  return pos_;
}

// Method that returns the velocity vector of the particle
arma::vec Particle::velocity(){
  return vel_;
}

*/
