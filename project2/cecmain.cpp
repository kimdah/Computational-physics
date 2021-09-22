
#include <iostream>
#include <armadillo>
#include <string>

#define pi 3.14159265359

using namespace std;
using namespace arma;


double find_max_value(arma::mat A, int& k, int& l);
void task_4b();
double find_max_value(); // ta bort?

arma::mat create_tridiagonal(const arma::vec& a, const arma::vec& d, const arma::vec& e);
arma::mat create_tridiagonal(int n, double a, double d, double e);
arma::mat create_symmetric_tridiagonal(int n, double a, double d);

arma::vec analytical_eigenvalues(arma::mat A);
arma::mat analytical_eigenvectors(arma::mat A);



void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l); // fra code snippets
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        const int maxiter, int& iterations, bool& converged);

void write_to_file(arma::mat output, String filename);

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


  task_4b(); //Solution to task 4b

  // 7 a - in separate function?
  int n = 10 // steps
  int N = n-1; // matrix size
  int h = 1./n; // stepsize

  // setting up scaled xhat:
  arma::vec xhat = arma::vec(n+1);
  xhat(0) = 0;
  for (int i = 1; i < n; i++){
    xhat(i) = xhat(i-1) + i*h;
  }
  xhat(n+1) = 1;


  arma::mat A = create_symmetric_tridiagonal(N, -1/(h*h), 2/(h*h));
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
// find three smallest values
arma::vec find_3min(arma::vec eigenvalues){
  // HOW TO DO THIS IN A SMART WAY?

  for (int i = 0; )
}

void write_to_file(arma::mat output, String filename){
  int width = 30;
  int prec = 10;
  ofstream ofile;
  ofile.open("output.txt"); //string(filename)
        << setw(width) << setprecision(prec) << scientific << 2 << endl; // temp
  ofile.close(); //close file
}






// Create tridiagonal matrix from vectors.
// - lower diagonal: vector a, lenght n-1
// - main diagonal:  vector d, lenght n
// - upper diagonal: vector e, lenght n-1
arma::mat create_tridiagonal(const arma::vec& a, const arma::vec& d, const arma::vec& e)
{
  int n = d.size();
  // Start from identity matrix
  arma::mat A = arma::mat(n, n, fill::eye);
  A = A*d; //not work?

  // Fill first row (row index 0)
  // OBS! ANTAR A(ROW, COLUMN)
  A(0,0) = d(0);
  A(0,1) = e(0);

  // Loop that fills rows 2 to n-1 (row indices 1 to n-2)
  for (int r = 1; r <= n-2; r++){
    A(r, r-1) = a(r-1);
    A(r, r) = d(r); // diagonal element
    A(r, r+1) = e(r);
  }

  // Fill last row (row index n-1)
  A(n-1, n-2) = a(n-2);
  A(n-1, n-1) = d(n-1);

  return A;
}


// Create a tridiagonal matrix tridiag(a,d,e) of size n*n
// from scalar input a, d and e
arma::mat create_tridiagonal(int n, double a, double d, double e)
{
  arma::vec a_vec = arma::vec(n-1, arma::fill::ones) * a;
  arma::vec d_vec = arma::vec(n, arma::fill::ones) * d;
  arma::vec e_vec = arma::vec(n-1, arma::fill::ones) * e;

  // Call the vector version of this function and return the result
  return create_tridiagonal(a_vec, d_vec, e_vec);
}


// Create a symmetric tridiagonal matrix tridiag(a,d,a) of size n*n
// from scalar input a and d.
arma::mat create_symmetric_tridiagonal(int n, double a, double d)
{
  // Call create_tridiagonal and return the result
  return create_tridiagonal(n, a, d, a);
}





arma::mat analytical_eigenvectors(arma::mat A){ // 3, vurder aa samle disse i 1 funk
  // Denne gir riktige verdier, men fortegnene er feil!!!
  int N = arma::size(A)(0);
  double d = A(0,0);
  double a = A(0,1);

  arma::mat v(N,N);

  for (int i = 1; i <= N; i++){
    for (int j = 1; j <= N; j++){
      v(j-1,i-1) = sin((i*j*pi)/(N+1)); // Aji fordi det i definisjonen var i som ga kolonnevektorene.
    }
  }
  return arma::normalise(v);

}
arma::vec analytical_eigenvalues(arma::mat A){ // 3
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

  //for loop runs through the non-diagonal matrix elements under the diagonal.
  for (int j=0; j<=N-1; j++){

          for (int i=1+j; i<=N-1; i++){
            if(abs(A(i,j))>abs(max_value)){
              max_value= A(i,j);
              k = j, l=i; //column k and row l

            }
      }
  }

  return max_value;
}

//---------------Task 4B-------------
void task_4b(){

  //defines and fills the matix shown in task 4b)
  arma::mat B_4 = arma::mat(4, 4).fill(0.);
  for (int i = 0; i < 3; i++){  // row
      for (int j = 0; j < 3; j++){ // // column
        if (i == j){
          B_4(i,j) = 1;
        }
      }
    }
    B_4(0,3) = 0.5; B_4(1,2) = -0.7; B_4(2,1) = -0.7; B_4(3,0) = 0.5;


  //prints the matrix to terminal and calls the function from task 3
  //to find en print the max value of non-diagonal matrix element in the
  //sub triangular matrix.

  int k; int l;
  cout<< B_4 <<endl;
  //returns max value and assigns k as the column index and l as the row index
  cout <<"max value: "<< find_max_value(B_4,k,l) <<" row: "<<l<<" column: "<< k << endl;
}

//---------------Task 4B(end)-------------
