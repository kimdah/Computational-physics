
#include <iostream>
#include <armadillo>

#define pi 3.14159265359

using namespace std;


double find_max_value();
arma::vec analytical_eigenvalues(arma::mat A);
arma::mat analytical_eigenvectors(arma::mat A);

int main(int argc, char const *argv[]) {

  // Problem 3
  int N = 6;    // size of matrix
  int n = N+1; // number of steps
  double h = 1./n; // stepsize
  arma::mat A = arma::mat(N, N).fill(0.);

  cout << "h = " << h << endl;

  for (int i = 0; i < N; i++){  // row
    for (int j = 0; j < N; j++){ // // column
      if (i == j){
        A(i,j) = 2./(h*h);
      } else if ((j == i+1) || (i == j+1)){
        A(i,j) = -1./ (h*h);
      }
    }
  }
  //cout << A;
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, A);

  cout << "Eigenvalues:\n" << eigval << endl; // printing out, delete later
  cout << "Eigenvectors:\n" << eigvec << endl;

  arma::vec eigval_analytic = analytical_eigenvalues(A);
  arma::mat eigvec_analytic = analytical_eigenvectors(A);

  cout << "Analytical eigenvalues:\n" << eigval_analytic << endl;

  cout << "Analytical eigenvectors:\n" << eigvec_analytic << endl;


  return 0;
}

arma::mat analytical_eigenvectors(arma::mat A){ // vurder aa samle disse i en
  int N = arma::size(A)(0);
  double d = A(0,0);
  double a = A(0,1);

  arma::mat v(N,N);

  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      v(i,j) = sin(((i+1)*pi)/(N+1)); // i eller j inni?
    }
  }
  return v;

}
arma::vec analytical_eigenvalues(arma::mat A){
  // A is tridiagonal (a,d,a)
  int N = arma::size(A)(0);
  double d = A(0,0);
  double a = A(0,1);

  arma::vec lambda(N);

  for (int i = 1; i <= N; i++){
    lambda(i-1) = d + 2*a*cos((i*pi)/(N+1));
  }
  return lambda;
}


double find_max_value(arma::mat A, int& k, int& l){

  double max_value = 0;
  int N = arma::size(A)(0); //i is dimension N of matrix A


  for (int i=0; i<=N+1; i++){
          for (int j=0; j<=N+1; j++){

            if(abs(A(i,j))>abs(max_value)){
              max_value= A(i,j);
              k = i, l=j;
            }
      }
      }

  return max_value;
}
