#ifndef CRANK_HPP
#define CRANK_HPP
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>

#include <armadillo>

#include <complex>
#include <cmath>

using namespace std::complex_literals; // to use imaginary number i | DEMANDS c++14!
using namespace std;
using namespace arma;

class Crank {
    public:
    sp_cx_mat A_, B_;
    cx_mat U_, U_empty;
    mat V_;
    cx_vec u_;
    double deltat_;
    bool poutput_; // Sets whether we want the probability output or real+imaginary
    int M_, t_; // size of total matrix
    complex<double> r_;

    Crank(double h, double deltat, double T, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double v_0, int slits);
    // Performs simulations based on parameter inputs

    int get_k_index(int i, int j, int M);
    cx_vec construct_u_vec(cx_mat U, bool normalise);
    mat make_potential_box(double v0);
    mat make_potential_single_slit(double v0);
    mat make_potential_double_slit(double v0);
    mat make_potential_triple_slit(double v0);
    void make_matrices(int M, double h, double deltat, mat V, complex<double> r);
    sp_cx_mat make_matrix(complex<double> r, cx_vec d);
    cx_vec time_step(cx_vec u);
    cx_mat make_insert_wavepacket(int M, double h, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y);
    void print();
    int to_file(string s);
    int fritjofs(string s);
    cx_cube run_simulation(int t);
    cx_mat col_to_mat(cx_vec u);
    double sum_probabilies(cx_mat U); 
    void output_probabilities(cx_cube R, string filename);

};
#endif
