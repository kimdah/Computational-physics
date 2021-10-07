#include <armadillo>
#include <iostream>


int main(int argc, char const *argv[]) {
  /* code */
  arma::vec a = arma::vec(2).fill(2);
  arma::vec b = arma::vec(2).fill(1);


  //std::cout << a+b;
  std::cout << "norm a: " << arma::norm(a);
  std::cout << "norm b: " << arma::norm(b);
  std::cout << "norm a-b: " << arma::norm(a-b);
  //std::cout << "\n abs: " << abs(a-b);

  // Evolve the system one time step (dt) using Euler-Cromer
  void PenningTrap::evolve_Euler_Cromer(double dt){
    for (int p = 0; p < particles_.size(); p++){
      arma::vec r = particle_[p].pos_;
      arma::vec v = particle_[p].vel_;
      double m = particle_[p].m_;

      v = v + dt * total_force(p)/m;
      x = x + dt*v;

    }

  }

  return 0;
}
