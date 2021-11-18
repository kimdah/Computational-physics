#include "../include/Ising.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

Ising::Ising(int lattice_side_length, double T, int seed, int ordered_spin) {
    L_ = lattice_side_length;
    N_ = L_*L_;
    T_ = T;

    //Initalise randomness with Mersenne Twister 19937 random number generator
    generator.seed(seed);
    // proposal_pdf_ = normal_distribution(0.0, 1.0);
    // lattice_uniform_distribution_ = uniform_int_distribution(0, L_-1);
    // up_or_down_spin_ = uniform_int_distribution(0, 1);

    normal_distribution<double> proposal_pdf_(0.0, 1.0);
    uniform_int_distribution<int> lattice_uniform_distribution_(0, L_-1);
    uniform_int_distribution<int> up_or_down_spin_(0, 1);


    std::vector<double> boltzmann_factors_ = calc_boltzmann_factors(T);
    if (ordered_spin != 0){
      generate_ordered_lattice(ordered_spin);
    } else {
      generate_unordered_lattice();
    }
}

void Ising::generate_ordered_lattice(int spin) {
    vector<vector<int>> lattice(L_, vector<int>(L_, spin));
    s_ = lattice;
}

void Ising::generate_unordered_lattice() {
     vector<vector<int>> lattice(L_, vector<int>(L_, 1));
     for (int i=0; i<L_; i++){
        for (int j=0; j<L_; j++){
            lattice[i][j] = lattice[i][j] - 2 * up_or_down_spin_(generator); // Generates a 1 or a -1
        }
    }
    s_ = lattice;
}

std::vector<std::vector<int>> Ising::run_metropolis_MCMC(){
  // running one MC cycle for sampling

  int epsilon = 0;
  int E = 0;
  double total_energy_per_cycle = 0;
  double total_magnetization_per_cycle = 0;
  double energy_of_new_state = calc_tot_energy_of_state(s_);
  for (int c = 0; c < N_; c++){ // one MC cycle; attempt N spin flips
    // flip random spin
    int randRow = rand() % L_;
    int randCol = rand() % L_;
    s_[randRow][randCol] *= -1;

    // examining surrounding spins to figure out index in boltzmann_factor vector
    // for computing the probability ratio
    int deltaE = s_[randRow][randCol] * 2 * (
                 s_[randRow][(randCol - 1 + L_) % L_] // Neighbour to the left
               + s_[randRow][(randCol + 1) % L_] // Neighbour to the right
               + s_[(randRow + 1) % L_][randCol] // Neighbour above
               + s_[(randRow - 1 + L_) % L_][randCol]); // Neighbour below

    // finding the index to use in Boltzmann
    int index;
    // boltzmann factor depends on flipping a +1 to -1, so the value will have
    // reverse index when a negative spin is flipped to positive.
    if(s_[randRow][randCol] == 1){ // if spin has been flipped to positive
      index = 5 - deltaE/4 + 2; // reversing index
    }else{
      index = deltaE/4 + 2;
    }

    // Acceptance ratio
    double probability_ratio = boltzmann_factors_[index];
    uniform_real_distribution<double> uniform_real(0.0, 1.0);
    double r = uniform_real(generator);
    //std::uniform_real_distribution<double> distribution(0.0,1.0);
    //double r = distribution()
    if (deltaE < 0 ){

    }
    else if(r > probability_ratio){ //Rejected spin-flip
      s_[randRow][randCol] *= -1; // flip back
    }
    else {
      // Accept spin configuration candidate
      //int totalenergy = totalenergy + deltaE; //
      //epsilon += totalenergy/N;
      //double magnetization += calc_tot_magnetization_of_state(s_current);
      //E += totalenergy; // metroplis
      energy_of_new_state += deltaE;

      total_magnetization_per_cycle += calc_tot_magnetization_of_state(s_);
    }
  }
  E = total_energy_per_cycle;
  total_energy_per_cycle = energy_of_new_state;

  //
  //
  // exp_val_eps_per_cycle = epsilon / acceptedstates; // <eps>
  // exp_val_eps_per_cycle_squared = pow(epsilon,2) / acceptedstates; //<eps^2>
  // exp_val_m_per_cycle = magnetization/acceptedstates; //<m>
  // exp_val_m_per_cycle_squared = pow(magnetization,2) / acceptedstates; //<m^2>
  // exp_val_E_per_cycle = E / acceptedstates; //<E>
  // exp_val_E_per_cycle_squared = pow(E,2) / acceptedstates;
  // heatcapacity_per_cycle = (1./acceptedstates)*(1./pow(T,2))*(exp_val_E_per_cycle_squared - pow(exp_val_E_per_cycle
  // ,2)); //C_v = 1/N_ 1/kbT^2 (<E^2>-<E>^2)
  return s_;
}



int Ising::calc_tot_energy_of_state(std::vector<std::vector<int> > s){
  // finding the energy of a particular spin configuration s
  int energy;
  for(int i=1 ; i<L_+1 ; i++){ //the first row will be the Lth row
    for(int j=1 ; j<L_+1 ; j++){ //the first column will be the Lth column
      int i_index = (i + L_)%L_;
      int j_index = (j + L_)%L_;
      energy += s[i_index][j_index]*s[i_index-1][j_index] + s[i_index][j_index]*s[i_index][j_index-1];
    }
  }
  return energy;
}

int Ising::calc_tot_magnetization_of_state(std::vector<std::vector<int>> s){
  int magnetization;
  for(int i=1 ; i<L_+1 ; i++){ //the first row will be the Lth row
    for(int j=1 ; j<L_+1 ; j++){ //the first column will be the Lth column
      magnetization +=s[i][j];
    }
  }
  return abs(magnetization);
}


std::vector<double> Ising::calc_boltzmann_factors(double T){
  double kB = 1.38064852 * pow(10, -23);
  double beta = 1. / (kB*T);
  vector<double> boltzmann_values;
  boltzmann_values.push_back(exp(-beta*(-8))); // 0 +1 spins
  boltzmann_values.push_back(exp(-beta*(-4))); // 1 +1 spins
  boltzmann_values.push_back(exp(-beta*(0))); // = 1? // 2 +1 spins
  boltzmann_values.push_back(exp(-beta*(4))); // 3 +1 spins
  boltzmann_values.push_back(exp(-beta*(8))); // 4 +1 spins
  return boltzmann_values;
}

void Ising::analytical_2x2(double T){
  double kB = 1.38064852 * pow(10, -23);
  double beta = 1. / (kB*T);

  double Z = 2*exp(beta*8) + 2*exp(-beta*8) + 12;
  double exp_val_epsilon = (4./Z) * (exp(-beta*8) - exp(beta*8));
  double exp_val_abs_mag = (2./Z) * (exp(beta*8) + 2);
  double heat_capacity = (32./(kB*pow(T,2)*Z))*(exp(-beta*8)-exp(beta*8)- (2./Z)*(exp(-beta*16)+exp(beta*16)-2));
  double susceptibility = (2/(kB*T*Z)) * (exp(8*beta) +1 - ((2./Z)*(exp(16*beta)+ 4*exp(8*beta)+4)));

}
