
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
int get_k_index(int i, int j, int M);
void make_matrices(int M, double h, double deltat, sp_cx_mat V, double r);
sp_cx_mat make_matrix(double r, cx_vec d);
cx_vec time_step(sp_cx_mat A, sp_cx_mat B, cx_vec u);
sp_cx_mat make_wavepacket(int M, double h, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y);
sp_cx_mat A; // glboal variables
sp_cx_mat B;


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
  //cout << "V:\n"<< V<< endl;

  //makes matrces A and B
  make_matrices(5, 0.1, 0.1, V, 2);
  sp_cx_mat U = make_wavepacket(5, 0.1, 0.1, 0.1, 0.2, 0.2, 0.1, 0.1);

  for (int i = 0; i<5; i++) {
    for (int j = 0; j<5; j++) {
      cout << setw(25) << U(i,j);
    }
    cout << endl;
  }
  //cout << U << endl;



}

//crates matrix U at time t=0
sp_cx_mat make_wavepacket(int M, double h, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y){

  sp_cx_mat U = sp_cx_mat(M, M); //Creates the matrix U

  //calculates non-boundary condtions
  for(int i =1; i< M-1; i++){
    for(int j =1;j< M-1; j++){
      double x = i*h;
      double y = j*h;
      U(i,j)= exp(-(pow(x-x_c,2)/(2*pow(sigma_x,2)))-(pow(y-y_c,2)/(2*pow(sigma_y,2))) + 1i*p_x*(x-x_c)+ 1i*p_y*(y-y_c));
    }
  }

  cx_double bc= (0,1); //boundary condition(Need to find correct) only works with imaginary != 0

  //cout << M <<endl;
  //Filling in boundary conditions
  for(int i=0; i < M+1; i+=(M-1)){ //Should change loop conditions to make more clear
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
// Translates matrix (i,j) index to column (k) index which have values from 1 to M-2
int get_k_index(int i, int j, int M){
//return ((i % (M - 1)) - 1) + (M - 2) * (j - 1);
  return ((j - 1) * (M - 2)) + (i - 1);

  }

// Task3
cx_vec time_step(sp_cx_mat A, sp_cx_mat B, cx_vec u){
  int m_size = B.size();
  // Try to get this to work for part 1
  //cx_vec b = affmul(B,u.t()); //Calculates Bu = b (maybe cross() instead?) (did not work).
  cx_vec b = cx_vec(m_size);

  // Try to optimise this
  //matrix multiplication Bu=b(instead of affmul()
  for(int i =0;i< m_size; i++){
    for(int j =0;j< m_size; j++){
      b(i) += (u(i).real()*B(i,j).real()-B(i,j).imag()*u(i).imag())+1i*(u(i).real()*B(i,j).imag()+u(i).imag()*B(i,j).real());
    }
  }
  // PArt 2
  return spsolve(A,b.t());	//spsolve assumes sparse matix. Solves the matrix eq. Ax=b, maybe solve() instead.
 }


// Makes specialized A and B matrices (Task 2.3)
void make_matrices(int M, double h, double deltat, sp_cx_mat V, double r){
  // assuming r is a real number
  int mat_size = pow(M-2,2);
  cx_vec a = cx_vec(mat_size);
  cx_vec b = cx_vec(mat_size);

  for(int k = 0 ; k < mat_size ; k++){
    double real = (deltat/2) * V(k,k).real(); // these work, but gives unnecessary work
    double img = (deltat/2) * V(k,k).imag();
    // a(k) = cx_double(1 + 4*r - img, real); // assuming r is real
    // b(k) = cx_double(1 - 4*r + img, -real);

    // We want these to work(behold til gruppetime paa torsdag):
    cout << "i*V: " << cx_double(V(k,k)) * 1.0i << endl;
    a(k) = (1 + 4*r + 1.0i*(deltat/2*cx_double(V(k,k))));
    b(k) = (1 - 4*r - 1.0i*(deltat/2*cx_double(V(k,k))));
  }

  A = make_matrix(-r, a); // made these global variables - maybe make a class instead?
  B = make_matrix(r,  b);

}

// Makes a matrix on the general form of A and B
sp_cx_mat make_matrix(double r, cx_vec d){
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
