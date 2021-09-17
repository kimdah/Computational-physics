#include <armadillo>

double find_max_value();


int main(int argc, char const *argv[]) {

  // Problem 3
  int N = 6;    // size of matrix
  arma::mat A = arma::mat(N, N).fill(0.);

  for (int i = 0; i < N; i++){
    
  }
  double A_ij = A(i,j); //Assign element (i,j) of the matrix A to A_ij.

  return 0;
}


double find_max_value(arma::mat A, int& k, int& l){

double max_value = 0;
N = arma::size(A)(0); //i is dimension N of matrix A


for (int i=0; i<=N+1; j++){
        for (int j=0; j<=N+1; j++){
          
          if(abs(A(i,j))>abs(max_value)){
            max_value= A(i,j)
            k = i, l=j;
          }
    }
    }

return max_value
}

