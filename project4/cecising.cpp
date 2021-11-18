#include "../include/Ising.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

Ising::Ising(int lattice_side_length, double T, int seed, bool generate_new_lattice) {
    make_new_lattice = generate_new_lattice;
    L = lattice_side_length;
    N = L*L;
    //Initalise randomness with Mersenne Twister 19937 random number generator
    generator.seed(seed);
    proposal_pdf = normal_distribution(0.0,1.0);
    lattice_uniform_distribution = uniform_int_distribution(0, L-1);
    up_or_down_spin = uniform_int_distribution(0, 1);

    boltzmann_factors = boltzmann_factor(T);

    if (make_new_lattice) { generate_unordered_lattice(); }

}

void Ising::generate_unordered_lattice() {
     vector<vector<int>> lattice(L, vector<int>(L, 1));
     for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            lattice[i][j] = lattice[i][j] - 2*up_or_down_spin(generator);
        }
    }
    s_current = lattice;
}

void Ising::generate_ordered_lattice(int spin) {
    vector<vector<int>> lattice(L, vector<int>(L, spin));

}

std::vector<std::vector<int>>Ising::run_metropolis_MCMC(){
  // running one MC cycle for sampling

  int epsilon = 0;
  int E = 0;
  double total_energy_per_cycle = 0;
  double total_magnetization_per_cycle = 0;
  double energy_of_new_state = calc_tot_energy_of_state();
  for (int c = 0; c < N; c++){ // one MC cycle; attempt N spin flips
    // flip random spin
    int randRow = rand() % L;
    int randCol = rand() % L;
    s_current[randRow][randCol] *= -1; // flip the spin

    // examining surrounding spins to figure out index in boltzmann_factor vector
    // for computing the probability ratio
    int deltaE = s_current[randRow][randCol] * 2 * (
                 s_current[randRow][(randCol - 1 + L) % L] // Neighbour to the left
               + s_current[randRow][(randCol + 1) % L]; // Neighbour to the right
               + s_current[(randRow + 1) % L][randCol] // Neighbour above
               + s_current[(randRow - 1 + L) % L][randCol]) // Neighbour below

    // finding the index to use in Boltzmann
    int index;
    // boltzmann factor depends on flipping a +1 to -1, so the value will have
    // reverse index when a negative spin is flipped to positive.
    if(s_current[randRow][randCol] == 1){ // if spin has been flipped to positive
      index = 5 - deltaE/4 + 2; // reversing index
    }else{
      index = deltaE/4 + 2;
    }

    // Acceptance ratio
    double probability_ratio = boltzmann_factors[index];
    double r = uniform_real_distribution(0.0, 1.0);
    if (deltaE < 0 ){

    }
    if(r > probability_ratio){ //Rejected spin-flip
      s_current[randRow][randCol] *= -1; // flip back
    }
    else{
      // Accept spin configuration candidate
      //int totalenergy = totalenergy + deltaE; //
      //epsilon += totalenergy/N;
      //double magnetization += calc_tot_magnetization_of_state(s_current);
      //E += totalenergy; // metroplis
      energy_of_new_state += deltaE;

      total_magnetization_per_cycle += calc_tot_magnetization_of_state(s_current);
    }
  }
  E = total_energy_per_cycle;
  total_energy_per_cycle = energy_of_new_state;



  exp_val_eps_per_cycle = epsilon / acceptedstates; // <eps>
  exp_val_eps_per_cycle_squared = pow(epsilon,2) / acceptedstates; //<eps^2>
  exp_val_m_per_cycle = magnetization/acceptedstates; //<m>
  exp_val_m_per_cycle_squared = pow(magnetization,2) / acceptedstates; //<m^2>
  exp_val_E_per_cycle = E / acceptedstates; //<E>
  exp_val_E_per_cycle_squared = pow(E,2) / acceptedstates;
  heatcapacity_per_cycle = (1./acceptedstates)*(1./pow(T,2))*(exp_val_E_per_cycle_squared - pow(exp_val_E_per_cycle
  ,2)); //C_v = 1/N 1/kbT^2 (<E^2>-<E>^2)
  return s_current;
}



int Ising::calc_tot_energy_of_state(std::vector<std::vector<int> > s){
  // finding the energy of a particular spin configuration s
  int energy;
  for(int i=1 ; i<L+1 ; i++){ //the first row will be the Lth row
    for(int j=1 ; j<L+1 ; j++){ //the first column will be the Lth column
      int i_index = (i + L)%L;
      int j_index = (j + L)%L;
      energy += s[i_index][j_index]*s[i_index-1][j_index] + s[i_index][j_index]*s[i_index][j_index-1];
    }
  }
  return energy;
}

int Ising::calc_tot_magnetization_of_state(std::vector<std::vector<int> > s){
  int magnetization;
  for(int i=1 ; i<L+1 ; i++){ //the first row will be the Lth row
    for(int j=1 ; j<L+1 ; j++){ //the first column will be the Lth column
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


  std::unordered_map<int, double> boltzmanns;
  boltzmanns[-8] = exp(-beta*(-8));
  return boltzmann_values;
}

void analytical_2x2(double T){
  double kB = 1.38064852 * pow(10, -23);
  double beta = 1. / (kB*T);

  double Z = 2*exp(beta*8) + 2*exp(-beta*8) + 12;
  double exp_val_epsilon = (4./Z) * (exp(-beta*8) - exp(beta*8));
  double exp_val_abs_mag = (2./Z) * (exp(beta*8) + 2)
  double heat_capacity = (32./(kB*pow(temp,2)*Z))*(exp(-beta*8)-exp(beta*8)- (2./Z)*(exp(-beta*16)+exp(beta*16)-2));
  double susceptibility =

}
