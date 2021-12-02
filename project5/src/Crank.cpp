#include "../include/Crank.hpp"
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

Crank::Crank(int M, double h, double deltat, double r, double v0) {
  M_ = M;

  V_ = make_potential(v0); // initialise V

  //makes matrices A and B
  make_matrices(M_, h, deltat, V_, r); // random variables!! change

  // Commented out to test errors:
  U_ = make_wavepacket(M_, h, 0.1, 0.1, 0.2, 0.2, 0.1, 0.1);
  u_ = construct_u_vec(U_,true);
  time_step(A_,B_,u_);

}

// Initialise potential V (time-independent) - make complex if necessary
mat Crank::make_potential(double v0){
  mat V = mat(M_,M_); //
  double infty = numeric_limits<double>::max(); // BETTER IDEA?
  V.col(0) = vec(M_, fill::value(infty));
  V.col(M_-1) = vec(M_).fill(infty);
  V.row(0) = rowvec(M_).fill(infty);
  V.row(M_-1) = rowvec(M_).fill(infty);
  //submat(first_row, first_col, last_row, last_col)
  V.submat(1, 1, M_-2, M_-2) = mat(M_-2,M_-2).fill(v0); // filling inner matrix
  return V;
}


sp_cx_mat Crank::make_wavepacket(int M, double h, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y){

  sp_cx_mat U = sp_cx_mat(M, M); //Creates the matrix U

  //calculates non-boundary condtions
  for(int i =1; i< M-1; i++){
    for(int j =1;j< M-1; j++){
      double x = i*h;
      double y = j*h;
      U(i,j)= exp(-(pow(x-x_c,2)/(2*pow(sigma_x,2)))-(pow(y-y_c,2)/(2*pow(sigma_y,2))) + 1i*p_x*(x-x_c)+ 1i*p_y*(y-y_c));
    }
  }

  cx_double bc= cx_double(0,0); //boundary condition(Need to find correct) only works with imaginary != 0

  //cout << M <<endl;
  //Filling in boundary conditions
  for(int i=0; i < M+1; i+=(M-1)){ // Should change loop conditions to make more clear
    for(int j=0; j< M; j++){
        U(i,j) = bc;
        U(j,i) = bc;
    }
  }

  /*
  Currently not working. Need to normalize U according to problem 4


  //finding normalisation constant
  cx_double normalization_cosntant = (0,0);
  for(int i=0;i<M; i++){
    for(int j=0;j<M; j++){
      normalization_cosntant += pow(abs(U(i,j)),2); //=p(x,y;t)?
    }
  }

  cout << normalization_cosntant << endl;
  //normalising
  for(int i=0;i<M; i++){
    for(int j=0;j<M; j++){
      U(i,j) = U(i,j)/normalization_cosntant;
    }
  }
  */
  return U;
}

// Problem 2-1
// Translates matrix (i,j) index to column (k) index which have values from 1 to M-1
int Crank::get_k_index(int i, int j, int M){
//return ((i % (M - 1)) - 1) + (M - 2) * (j - 1);
  return ((j - 1) * (M - 2)) + (i - 1);

}

// constructs the u vector based on U matrix
cx_vec Crank::construct_u_vec(sp_cx_mat U, bool normalise){
  int M = sqrt(U.size()); //size() gives M², assumes U to be quadratic
  cx_vec u = cx_vec(pow(M-2,2));

  for(int i = 1; i< M-1; i++){
    for(int j = 1; j< M-1; j++){
      u(get_k_index(i, j, M)) = U(i, j);
    }
  }
  if(normalise){ //Normalizes u, but not sure if correct
    sp_cx_mat normalisation_factor;
    normalisation_factor = u.t()*u; //this is a 1x1 matrix
    u = u/normalisation_factor(0, 0);
  }
  return u;
}

// Task3
cx_vec Crank::time_step(sp_cx_mat A, sp_cx_mat B, cx_vec u){
  cx_vec a = B * u;
  return spsolve(A,a);

  // int m_size = sqrt(B.size());  //assuming quadratic matrix
  // // Try to get this to work for part 1
  // //cx_vec b = affmul(B,u.t()); //Calculates Bu = b (maybe cross() instead?) (did not work).
  // cx_vec b = cx_vec(m_size);

  // // Try to optimise this
  // //matrix multiplication Bu=b(instead of affmul()
  // for(int i =0;i< m_size; i++){
  //   for(int j =0;j< m_size; j++){
  //     b(i) += (u(i).real()*B(i,j).real()-B(i,j).imag()*u(i).imag())+1i*(u(i).real()*B(i,j).imag()+u(i).imag()*B(i,j).real());
  //   }
  // }

  // PArt 2

  //return spsolve(A,b);	//spsolve assumes sparse matix. (removed .t())
 }


// Makes specialized A and B matrices (Task 2.3)
void Crank::make_matrices(int M, double h, double deltat, mat V, double r){
  // assuming r is a real number
  int mat_size = pow(M-2,2);
  cx_vec a = cx_vec(mat_size);
  cx_vec b = cx_vec(mat_size);

  // Wrong indexing for V of size M
  // for(int k = 0 ; k < mat_size ; k++){
  //   //double real = (deltat/2) * V(k,k).real(); // these work, but gives unnecessary work
  //   //double img = (deltat/2) * V(k,k).imag();
  //   // a(k) = cx_double(1 + 4*r - img, real); // assuming r is real
  //   // b(k) = cx_double(1 - 4*r + img, -real);
  //
  //   // We want these to work(behold til gruppetime paa torsdag):
  //   //cout << "i*V: " << cx_double(V(k,k)) * 1.0i << endl;
  //   a(k) = (1 + 4*r + 1.0i*(deltat/2*cx_double(V(k,k))));
  //   b(k) = (1 - 4*r - 1.0i*(deltat/2*cx_double(V(k,k))));
  // }

  // Alternative:
  for (int i = 1; i < M-1; i++){ // Excluding boundaries in V (infinity) - is that ok?
    for (int j = 1; j < M-1; j++){
      int k = get_k_index(i,j,M);
      a(k) = (1 + 4*r + 1.0i*(deltat/2*cx_double(V(i,j))));
      b(k) = (1 - 4*r - 1.0i*(deltat/2*cx_double(V(i,j))));
    }
  }

  A_ = make_matrix(-r, a); // made these global variables - maybe make a class instead?
  B_ = make_matrix(r,  b);
}

// Makes a matrix on the general form of A and B
sp_cx_mat Crank::make_matrix(double r, cx_vec d){
  int S = d.size(); // (M-2)^2
  int s = sqrt(S); // (M-2)
  sp_cx_mat M = sp_cx_mat(S, S);

  // Making diag
  for (int i = 0; i < S; i+=s){ // last index = S-s
    sp_cx_mat D(s, s);
    D.diag() = d.subvec(i, i+s-1);// diagonal
    D.diag(-1) = cx_vec(s-1).fill(r); // subdiagonal
    D.diag(1) = cx_vec(s-1).fill(r); // superdiagonal
    //submat(first_row, first_col, last_row, last_col)
    M.submat(i, i, i+s-1, i+s-1) = D; //ex: (0,0,2,2), (3,3,5,5)
  }

  // Making non-diag, non-corners
  for (int i = s; i < S; i+=s){ // last index= S-s
    sp_cx_mat ND(s,s);
    ND.diag() = cx_vec(s).fill(r); // fill diagonal with r value
    //submat(first_row, first_col, last_row, last_col)
    M.submat(i-s, i, i-1, i+s-1) = ND; //ex: i=s=3: (0,3,2,5) i=2s=6: (3,6,5,8)
    M.submat(i, i-s, i+s-1, i-1) = ND; //ex: i=s=3: (3,0,5,2) i=2s=6: (6,3,8,5)
  }

  // Corners are 0 matrix, which they are already initialized as through sp_ (sparse)
  return M;
}

void Crank::print() {
    for (int i = 0; i<5; i++) {
        for (int j = 0; j<5; j++) {
        cout << setw(25) << U_(i,j);
        }
        cout << endl;
    }
}
