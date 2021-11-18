#ifndef ISING_HPP
#define ISING_HPP
#include <random>

class Ising {
    public:
        int L, N, seed;
        double T, totalenergy, exp_val_eps_per_cycle, exp_val_eps_per_cycle_squared,
        exp_val_m_per_cycle, exp_val_m_per_cycle_squared, exp_val_E_per_cycle,
        exp_val_E_per_cycle_squared, heatcapacity_per_cycle;
        bool make_new_lattice;
        mt19937 generator;
        std::vector<double> calc_boltzmann_factors;
        std::vector<std::vector<int>> s_current;
        normal_distribution<double> proposal_pdf;
        uniform_int_distribution<int> lattice_uniform_distribution, up_or_down_spin;
        void generate_unordered_lattice();
        void generate_ordered_lattice(int spin);


    Ising(int lattice_side, double T, int seed);


}











#endif
