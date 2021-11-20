#include "../include/Ising.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iomanip>

using namespace std;

Ising::Ising(int lattice_side_length, double T, int seed, int ordered_spin) {
    L_ = lattice_side_length;
    N_ = L_*L_;
    T_ = T;
    mag_per_spin_ = 0;
    epsilon_ = 0;
    tot_cycles_ = 0;

    //Initalise randomness with Mersenne Twister 19937 random number generator
    generator_.seed(seed);
    proposal_pdf_ = normal_distribution<double>(0.0, 1.0);
    uniform_real_ = uniform_real_distribution<double>(0.0, 1.0);
    lattice_uniform_distribution_ = uniform_int_distribution<int>(0, L_-1);
    up_or_down_spin_ = uniform_int_distribution<int>(0, 1);

    boltzmann_factors_ = calc_boltzmann_factors(T);

    if (ordered_spin != 0){
      generate_ordered_lattice(ordered_spin);
    } else {
      generate_unordered_lattice();
    }
    calc_tot_magnetization_of_state();
    calc_energy_of_lattice_state();
}

void Ising::generate_ordered_lattice(int spin) {
    vector<vector<int>> lattice(L_, vector<int>(L_, spin));
    s_ = lattice;
}

void Ising::generate_unordered_lattice() {
     vector<vector<int>> lattice(L_, vector<int>(L_, 1));
     for (int i=0; i<L_; i++){
        for (int j=0; j<L_; j++){
            lattice[i][j] = lattice[i][j] - 2 * up_or_down_spin_(generator_); // Generates a 1 or a -1
        }
    }
    s_ = lattice;
}

vector<vector<int>> Ising::run_metropolis_MCMC(){
  int randRow, randCol, index, deltaE;

  for (int c = 0; c < N_; c++){ // one MC cycle; attempt N spin flips
    // flip random spin
    randRow = lattice_uniform_distribution_(generator_);
    randCol = lattice_uniform_distribution_(generator_);
    //s_[randRow][randCol] *= -1; // Flipping a random spin. Will get flipped back if it fails later

    // examining surrounding spins to figure out index in boltzmann_factor vector
    // for computing the probability ratio. Computes the difference in energy after flipping a spin
    deltaE = s_[randRow][randCol] * 2 // The spin flip happens here (no minus).
           *(s_[randRow][(randCol - 1 + L_) % L_] // Neighbour to the left
           + s_[randRow][(randCol + 1) % L_] // Neighbour to the right
           + s_[(randRow + 1) % L_][randCol] // Neighbour above
           + s_[(randRow - 1 + L_) % L_][randCol]); // Neighbour below

    // boltzmann factor depends on flipping a +1 to -1, so the value will have
    // reverse index when a negative spin is flipped to positive.

    // if(s_[randRow][randCol] == 1){ // if spin has been flipped to positive
    //   index = 5 - deltaE/4 + 2; // reversing index
    // }else{
    //   index = deltaE/4 + 2;
    // }
    index = deltaE/4 + 2;
    // Acceptance ratio
    double probability_ratio = boltzmann_factors_[index]; // w_i/w_j = exp(-beta*deltaE)

    // Use uniform or normal? Wasn't normal suggested in lectures?
    //double r = uniform_real_(generator_);
    double r = uniform_real_(generator_);

    if (r <= probability_ratio ){ //|| deltaE < 0
      // Accept spin configuration candidate
      // Always accept for energy reducing flips
      s_[randRow][randCol] *= -1;
      totalenergy_ += deltaE;
      magnetisation_ += 2 * s_[randRow][randCol]; // Equation 13.7 in lectures2015 M_(i+1) = M_i + 2*s_(i+1) (= +/- 2 )
    } else {
      //Reject spin-flip
      //s_[randRow][randCol] *= -1;
    }
  }
  //Adding the values from each cycle, so it can be used to find exp values.
  epsilon_ += totalenergy_/N_;
  mag_per_spin_ += 1.0*magnetisation_/ N_;
  accumulatedtotalenergy_ += totalenergy_;
  accumulatedtotalmagnetization_ += magnetisation_;
  tot_cycles_ += 1;
  return s_; // not neccessary to return s_?
}

double Ising::mean(double value, int n_cycles){
  return value / n_cycles;
}

double Ising::expval_epsilon(int n_cycles){
  // totalenergy er ikke per cycle, burde den vaere det?
  //return totalenergy / N_;
  return mean(epsilon_, n_cycles);
}

double Ising::expval_mag_per_spin(int n_cycles){
  //return magnetisation / N_;
  return mean(abs(mag_per_spin_), n_cycles);
}

double Ising::heat_capacity(int n_cycles){
  return (1./N_)*(1./n_cycles)*(1./pow(T_,2))*(mean(pow(accumulatedtotalenergy_, 2), n_cycles) - pow(mean(accumulatedtotalenergy_, n_cycles), 2)); //C_v = 1/N_ 1/kbT^2 (<E^2>-<E>^2)
}

double Ising::susceptibility(int n_cycles){
  return (1./N_)*(1./n_cycles)*(1./T_)*(mean(pow(accumulatedtotalmagnetization_, 2), n_cycles) - pow(mean(accumulatedtotalmagnetization_, n_cycles), 2));
}


// Working
void Ising::calc_energy_of_lattice_state() {
  double energy = 0;
  for (int i=0; i<L_; i++){
    for (int j=0; j<L_; j++){
      energy +=  - s_[i][j] * s_[(i+1)%L_][j]  +  s_[i][j] * s_[i][(j+1)%L_];
    }
  }
  totalenergy_ = energy;
}
// Overload of function above
// Working
int Ising::calc_energy_of_lattice_state(vector<vector<int> > s) {
  int energy = 0;
  for (int i=0; i<L_; i++){
    for (int j=0; j<L_; j++){
      energy +=  s[i][j] * s[(i+1)%L_][j]  +  s[i][j] * s[i][(j+1)%L_];
    }
  }
  totalenergy_ = energy;
  return energy;
}

/* // Not working. Gives segmentation fault
int Ising::calc_tot_energy_of_state(vector<vector<int> > s){
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
 */

// Working
void Ising::calc_tot_magnetization_of_state(){
  int magnetisation = 0;
  for(int i=0 ; i<L_ ; i++){ //the first row will be the Lth row
    for(int j=0 ; j<L_ ; j++){ //the first column will be the Lth column
      magnetisation +=s_[i][j];
    }
  }
  magnetisation_ = magnetisation;
}

int Ising::calc_tot_magnetization_of_state(vector<vector<int>> s){
  int magnetisation = 0;
  for(int i=0 ; i<L_ ; i++){ //the first row will be the Lth row
    for(int j=0 ; j<L_ ; j++){ //the first column will be the Lth column
      magnetisation += s_[i][j];
    }
  }
  return abs(magnetisation);
}

// Working. Things make sense
vector<double> Ising::calc_boltzmann_factors(double T){
  //double kB = 1.38064852 * pow(10, -23);
  //double beta = 1. / (kB*T);
  double beta = 1. / (T);
  vector<double> boltzmann_values;
  boltzmann_values.push_back(exp(-beta*(-8))); // 0 +1 spins
  boltzmann_values.push_back(exp(-beta*(-4))); // 1 +1 spins
  boltzmann_values.push_back(exp(-beta*(0))); // = 1? // 2 +1 spins
  boltzmann_values.push_back(exp(-beta*(4))); // 3 +1 spins
  boltzmann_values.push_back(exp(-beta*(8))); // 4 +1 spins
  return boltzmann_values;
}




// Adds current parameters to referenced ofstream
void Ising::write_parameters_to_file(ofstream& ofile) {
  int width = 16;
  int prec  = 8;

  ofile << setw(width) << tot_cycles_;
  ofile << setw(width) << totalenergy_;
  ofile << setw(width) << magnetisation_;
  ofile << setw(width) << expval_epsilon(tot_cycles_);
  ofile << setw(width) << expval_mag_per_spin(tot_cycles_);
  ofile << setw(width) << heat_capacity(tot_cycles_);
  ofile << setw(width) << susceptibility(tot_cycles_);
  ofile << endl;
  sample_+=1;
}

void Ising::print(){
  cout << "Energy of lattice: " << totalenergy_ << "\n";
  cout << "Magnetisation of lattice: " << magnetisation_ << "\n";
    for (int i=0; i<L_; i++){
        for (int j=0; j<L_; j++){
          if(s_[i][j] > 0) {
            cout << "+";
          }
            cout << s_[i][j] << " ";
        }
        cout << "\n";
    }
}
