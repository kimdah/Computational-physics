#ifndef ISING_HPP
#define ISING_HPP
#include <random>

class Ising {
    public:
        int L_, N_, seed;
        double T_, totalenergy_;
        mt19937 generator;
        std::vector<double> boltzmann_factors_;
        std::vector<std::vector<int>> s_;
        normal_distribution<double> proposal_pdf_;
        uniform_int_distribution<int> lattice_uniform_distribution_, up_or_down_spin_;

    Ising(int lattice_side, double T, int seed);


}











#endif