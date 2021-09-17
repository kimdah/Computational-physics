
#include <iostream>
#include <armadillo>

using namespace std;


double find_max_value();


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




  return 0;
}

double analytical_eigenvectors(){

}
double analytical_eigenvalues(){
  
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
