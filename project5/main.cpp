
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <complex>
#include <armadillo>


using namespace std;
using namespace arma;
using namespace std::complex_literals;
// Performs simulations based on parameter inputs
int change_index(int i, int j, int M);
//void make_matrices(double r, cx_vec a, cx_vec b);
void make_matrices(int M, double h, double deltat, sp_cx_mat V, double r);
sp_cx_mat make_matrix(double r, cx_vec d);
//rowvec time_step(rowvec u);



int main(int argc, char const *argv[]) {
  // int M=5;
  // double deltat;
  // double h;
  // sp_cx_mat A= sp_cx_mat((m-2)^2, (m-2)^2).fill:zeros;
  // sp_cx_mat B= sp_cx_mat((m-2)^2, (m-2)^2).fill:zeros;
  // sp_cx_mat V= sp_cx_mat((m-2)^2, (m-2)^2).fill:zeros;
  // cx_vec a = cx_vec((m-2)^2-1);
  // cx_vec b = cx_vec((m-2)^2-1);
  //
  // for(int k=0 ; k < (m-2)^2 ; k++){
  //   a(k) = (1 + 4*r, deltat/2 * V(k,k);
  //   b(k) = (1 - 4*r, -deltat/2 * V(k,k);
  //   A(k,k) = a(k);
  //   B(k,k) = b(k);
  // }

//cout << change_index(3,3,5)<<endl;

  // // Testing make_matrix methods (task 2.2)
  // cx_vec aa(9, fill::ones); // fill::randu
  cx_vec t(2*2, fill::value(3));
  //make_matrix(2, t);
  sp_cx_mat V = sp_cx_mat(9,9);
  V.diag() = cx_vec(9, fill::ones);


  make_matrices(5, 0.1, 0.1, V, 2);

}

//changes index for the u vector(column), i and j can have values from 0 to M-2
int change_index(int i, int j, int M){return ((i%(M-1))-1)+ (M-2)*(j-1);}

// rowvec time_step(rowvec u){
// 	rowvec b = affmul(B,u); //Calculates Bu = b (maybe cross() instead?)
// 	return spsolve(A,b);	//spsolve assumes sparse matix, maybe solve() instead.
// }


// Makes specialized A and B matrices (2.3)
void make_matrices(int M, double h, double deltat, sp_cx_mat V, double r){
  int mat_size = pow(M-2,2);
  cx_vec a = cx_vec(mat_size);
  cx_vec b = cx_vec(mat_size);

  for(int k = 0 ; k < mat_size ; k++){
    cout << "V real: " << V(k,k).real();
    double real = (deltat/2) * V(k,k).real();
    cx_double img = (deltat/2) * V(k,k).imag();
    cout << "real: "<< real << "\nimag: "<< img << endl;
    //cx_double img = (deltat/2) * V(k,k);
    //cx_double img = (deltat/2) * V(k,k);
    //cout << "V: "<< V(k,k) << "\nimg: " << endl;
    //a(k) = (1 + 4*r + 1i*((deltat/2)*V(k,k))); // (real, imaginary)
    //b(k) = (1 - 4*r + 1i*((-deltat/2)*V(k,k)));
  }

  sp_cx_mat A = make_matrix(-r, a);
  sp_cx_mat B = make_matrix(r,b);
  cout << "A:\n" << A << "\nB:\n" << B << endl;
  // perhaps make it return the matrices somehow?

}

// Makes a matrix on the general form of A and B
sp_cx_mat make_matrix(double r, cx_vec d){
  int S = d.size(); // (M-2)^2
  int s = sqrt(S); // (M-2)
  sp_cx_mat M = sp_cx_mat(S, S);

  // Making diag
  for (int i = 0; i < S; i+=s){ // last index i = S-s
    sp_cx_mat D(s,s);
    D.diag() = d.subvec(i,i+s-1);//cx_vec(s,fill::value(d(i))); // diagonal
    D.diag(-1) = cx_vec(s-1,fill::value(r)); // subdiagonal
    D.diag(1) = cx_vec(s-1,fill::value(r)); // superdiagonal
    M.submat(i,i,i+s-1,i+s-1) = D; //submat(first_row, first_col, last_row, last_col)
  }

  // Making non-diag, non-corners
  for (int i = s; i < S; i+=s){ // last index i = S-s
    sp_cx_mat ND(s,s);
    ND.diag() = cx_vec(s,fill::value(8)); // fill diagonal with r value
    //submat(first_row, first_col, last_row, last_col)
    M.submat(i-s,i,i-1,i+s-1) = ND; //feks: i=s=3: (0,3,2,5) i=2s=6: (3,6,5,8)
    M.submat(i, i-s, i+s-1, i-1) = ND; // feks i=s=3: (3,0,5,2) i=2s=6: (6,3,8,5)
  }

  // Corners are 0 matrix, which they are already initialized as through sp_ (sparse)

  return M;
}
