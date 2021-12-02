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
    sp_cx_mat A_, B_, U_;
    mat V_;
    cx_vec u_;
    int M_; // size of total matrix

    Crank(int M, double h, double deltat, double r, double v0);
    // Performs simulations based on parameter inputs

    int get_k_index(int i, int j, int M);
    cx_vec construct_u_vec(sp_cx_mat U, bool normalise);
    mat make_potential(double v0);
    void make_matrices(int M, double h, double deltat, mat V, double r);
    sp_cx_mat make_matrix(double r, cx_vec d);
    cx_vec time_step(sp_cx_mat A, sp_cx_mat B, cx_vec u);
    sp_cx_mat make_wavepacket(int M, double h, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y);
    void print();

};
#endif
