# include "Particle.hpp"




Particle::Particle(double q_in, double m_in, arma::vec pos_in, arma::vec vel_in){
  q_ = q;         // charge of particle, [e] elementary charge
  m_ = m_in;      // mass of particle in atomic mass unit u
  pos_ = pos_in;  // position vector r
  vel_ = vel_in;  // velocity vector v

}
