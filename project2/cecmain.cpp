
#include <iostream>
#include <armadillo>
#include <string>

#define pi 3.14159265359

using namespace std;
using namespace arma;


arma::mat create_tridiagonal(const arma::vec& a, const arma::vec& d, const arma::vec& e);
arma::mat create_tridiagonal(int N, double a, double d, double e);
arma::mat create_symmetric_tridiagonal(int N, double a, double d);






int main(int argc, char const *argv[]) {

  //problem 3
  int N = 6;                          //size of matrix NxN
  int n = N+1;                        //steps in matrix
  double a= -1./ ((1./n)*(1./n));     //super and sub diagonal elements
  double d= 2./((1./n)*(1./n));       //diagonal elements
  arma::mat A = create_symmetric_tridiagonal(N, a, d);


  // 7 a - in separate function?
  int n = 10; // steps
  int N = n-1; // matrix size
  double h = 1./n; // stepsize

  // setting up scaled xhat: // SAVE THIS!
  /*
  arma::vec xhat = arma::vec(n+1);
  xhat(0) = 0;
  for (int i = 1; i < n; i++){
    xhat(i) = xhat(i-1) + i*h;
  }
  xhat(n) = 1;
*/

  arma::mat A = create_symmetric_tridiagonal(N, -1./(h*h), 2./(h*h));
  cout << A;
  /*
  arma::vec eigenvalues = arma::vec(N);
  arma::mat eigenvectors = arma::mat(N,N);
  jacobi_eigensolver(A, double eps, eigenvalues, eigenvectors,
                        const int maxiter, int& iterations, bool& converged);

  arma::vec min3_index = find_3min(eigenvalues);
  for (int i = 0;  i < min3_index.size(); i++){
    int index = min3_index(i);
    arma::vec eigenvec = eigenvectors(:,i); // vil ha kolonnen her - sjekk at ok.

  }
  */



  return 0;
}

arma::mat create_symmetric_tridiagonal(int N, double a, double d) // SPECIAL VERSION
{

    // Problem 3
  int n = N+1;   //number of steps

  double h = 1./n; // stepsize
  arma::mat A = arma::mat(N, N).fill(0.);

  cout << "h = " << h << endl;

  for (int i = 0; i < N; i++){  // row
    for (int j = 0; j < N; j++){ // // column
      if (i == j){
        A(i,j) = d;
      } else if ((j == i+1) || (i == j+1)){
        A(i,j) = a;
      }
    }
  }
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, A);

  cout << "Eigenvalues:\n" << eigval << endl; // printing out, delete later
  cout << "Eigenvectors:\n" << eigvec << endl;

  arma::vec eigval_analytic = analytical_eigenvalues(A);
  arma::mat eigvec_analytic = analytical_eigenvectors(A);

  cout << "Analytical eigenvalues:\n" << eigval_analytic << endl;

  cout << "Analytical eigenvectors:\n" << eigvec_analytic << endl;

  return A;
}
/*
// find three smallest values
arma::vec find_3min(arma::vec eigenvalues){
  // HOW TO DO THIS IN A SMART WAY?

  for (int i = 0; )
}
*/

/*
void write_to_file(arma::mat output, String filename){
  int width = 30;
  int prec = 10;
  ofstream ofile;
  ofile.open("output.txt"); //string(filename)
        << setw(width) << setprecision(prec) << scientific << 2 << endl; // temp
  ofile.close(); //close file
}
*/
