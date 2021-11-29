
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


// Performs simulations based on parameter inputs
int change_index(int i, int j, int M);
void make_matrices(int M, double h, double deltat, sp_cx_mat V, double r);
sp_cx_mat make_matrix(double r, cx_vec d);
cx_vec time_step(sp_cx_mat A, sp_cx_mat B, cx_vec u);



int main(int argc, char const *argv[]) {
  // // Testing make_matrix methods (task 2.2)
  //cx_vec aa(9, fill::ones); // fill::randu
  // cx_vec t(2*2, fill::value(3));
  // make_matrix(2, t);


  // Some of code demands Armadillo version 10.6 - check your version here:
  // arma::arma_version ver;
	// std::cout << "ARMA version: "<< ver.as_string() << std::endl;

  // Testing make_matrices
  sp_cx_mat V = sp_cx_mat(9,9);
  // making it a bit random, but still sparse
  V.diag() = cx_vec(9, fill::randu);
  V.diag(-3) = cx_vec(9-3, fill::randu); // subdiagonal 3
  V.diag(2) = cx_vec(9-2, fill::randu); // superdiagonal 2
  cout << "V:\n"<< V<< endl;

  make_matrices(5, 0.1, 0.1, V, 2);

}

//changes index for the u vector(column), i and j can have values from 0 to M-2
int change_index(int i, int j, int M){return ((i%(M-1))-1)+ (M-2)*(j-1);}

// commented out to make things compile:
//Task3
 cx_vec time_step(sp_cx_mat A, sp_cx_mat B, cx_vec u){
 	int M = B.size();
 	//cx_vec b = affmul(B,u.t()); //Calculates Bu = b (maybe cross() instead?)
 	cx_vec b = cx_vec(pow(M-2,2));
 	//matrix multiplication Bu=b
 	for(int i =0;i< pow(M-2,2); i++){
 		for(int j =0;j< pow(M-2,2); j++){
 			b(i) += (u(i).real()*B(i,j).real()-B(i,j).imag()*u(i).imag())+1i*(u(i).real()*B(i,j).imag()+u(i).imag()*B(i,j).real());
 		}
 	}	
 	return spsolve(A,b.t());	//spsolve assumes sparse matix. Solves the matrix eq. Ax=b, maybe solve() instead.
 	}


// Makes specialized A and B matrices (2.3)
void make_matrices(int M, double h, double deltat, sp_cx_mat V, double r){
  // assuming r is real
  int mat_size = pow(M-2,2);
  cx_vec a = cx_vec(mat_size);
  cx_vec b = cx_vec(mat_size);

  for(int k = 0 ; k < mat_size ; k++){
    double real = (deltat/2) * V(k,k).real(); // trying these instead; these work, but gives unnecessary work
    double img = (deltat/2) * V(k,k).imag();
    //cout << "V(k,k)" << V(k,k)<<endl;
    //cout << "real: " << real << "     imag: " << img << endl;
    //cout << "i*V: " << i*V(k,k) << endl; // doesnt work for some reason!
    a(k) = cx_double(1 + 4*r - img, real); // assuming r is real
    b(k) = cx_double(1 - 4*r + img, -real);

    // We want these to work:
    //a(k) = (1 + 4*r + 1i*(deltat/2*V(k,k));
    //b(k) = (1 - 4*r - 1i*(deltat/2*V(k,k));
  }

  sp_cx_mat A = make_matrix(-r, a);
  sp_cx_mat B = make_matrix(r,b);
  //cout << "A:\n" << A << "\nB:\n" << B << endl;
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
    D.diag(-1) = cx_vec(s-1).fill(r); // subdiagonal
    D.diag(1) = cx_vec(s-1).fill(r); // superdiagonal
    //submat(first_row, first_col, last_row, last_col)
    M.submat(i,i,i+s-1,i+s-1) = D; //ex: (0,0,2,2), (3,3,5,5)
  }

  // Making non-diag, non-corners
  for (int i = s; i < S; i+=s){ // last index i = S-s
    sp_cx_mat ND(s,s);
    ND.diag() = cx_vec(s).fill(r); // fill diagonal with r value
    //submat(first_row, first_col, last_row, last_col)
    M.submat(i-s,i,i-1,i+s-1) = ND; //ex: i=s=3: (0,3,2,5) i=2s=6: (3,6,5,8)
    M.submat(i, i-s, i+s-1, i-1) = ND; //ex: i=s=3: (3,0,5,2) i=2s=6: (6,3,8,5)
  }

  // Corners are 0 matrix, which they are already initialized as through sp_ (sparse)
  return M;
}
