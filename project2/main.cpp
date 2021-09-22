
#include <iostream>
#include <armadillo>

#define pi 3.14159265359

using namespace std;


double find_max_value(arma::mat A, int& k, int& l);
void task_4b();

arma::mat create_symmetric_tridiagonal(int N, double a, double d);

double find_max_value();
arma::vec analytical_eigenvalues(arma::mat A);
arma::mat analytical_eigenvectors(arma::mat A);

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l, int N); // fra code snippets (why ref A?)
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        const int maxiter, int& iterations, bool& converged);


int main(int argc, char const *argv[]) {

  //problem 3
  int N = 6;                          //size of matrix NxN
  int n = N+1;                        //steps in matrix
  double a= -1./ ((1./n)*(1./n));     //super and sub diagonal elements
  double d= 2./((1./n)*(1./n));       //diagonal elements
  arma::mat A = create_symmetric_tridiagonal(N, a, d);

  task_4b(); //Solution to task 4b
  
  //Task 5:
  double epsilon = 1.0e-8; //Tolerance
  double max_number_iterations = (double) N * (double) N * (double) N;
  int iterations = 0;
  double max_value = find_max_value( A, &k, &l);
  arma::mat R = arma::mat( N, N, arma::fill::eye); //initializing R
    
  while ( fabs(max_value) > epsilon && (double) iterations < max_number_iterations ) {
      max:value = find_max_value( A, &k, &l);
      jacobi_rotate( A, R, k, l, N);
      iterations++;
  }
  cout << "Number of iterations: " << iterations << "\n";

  int number_of_rotations; //describes the number of rotations completed by jacobi_rotate()

  return 0;
}

arma::mat create_symmetric_tridiagonal(int N, double a, double d)
{

    // Problem 3
  int n = N+1;   //number of steps

  double h = 1./n; // stepsize
  arma::mat A = arma::mat(N, N).fill(0.);

  cout << "h = " << h << endl;

  for (int i = 0; i < N; i++){  // row
    for (int j = 0; j < N; j++){ // // column
      if (i == j){
        A(i,j) = 2./(h*h);
      } else if ((j == i+1) || (i == j+1)){
        A(i,j) = a;
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
    
  return A;
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

//---------------Task 5A------------------
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l, int N){
    //Computing tan (t), cos (c) and sin (s)
    double s, c;
    if ( A[k][l] != 0.0){
        double t, tau;
        tau = (A[l][l]-A[k][k])/(2*A[k][l]);
        if ( tau > 0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else {
            t = -1.0/( -tau + sqrt(1.0 + tau*tau));
        }
        c = 1.0/(sqrt(1+t*t));
        s = c*t;
    }
    else {
        c = 1.0;
        s = 0.0;
    }
    
    //Transform current A matrix
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A[k][k];
    a_ll = A[l][l];
    A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
    A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
    A[k][l] = 0.0;
    A[l][k] = 0.0;
    for ( int i = 0; i < N; i++ ) {
        if ( i != k && i != l ) {
            a_ik = A[i][k];
            a_il = A[i][l];
            A[i][k] = c*a_ik - s*a_il;
            A[k][i] = A[i][k];
            A[i][l] = c*a_il + s*a_ik;
            A[l][i] = A[i][l];
        }
        //Compute new eigenvectors
        r_ik = R[i][k];
        r_il = R[i][l];
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;
    }
    return;
}



//-----------Task 6-------------


void jacobi_scaling(int& number_of_rotations, int& a, int& d){

arma::mat A;
for (int N = 3; N < 100; N++){  
    A = create_symmetric_tridiagonal(N,a,d); //creates an NxN tridaiag symmetric matrix
    //jacobi_rotate(A);
    
    //should work as lons as Jacobi_rotate takes matrix A as an input
    //and updates a variable number_of_rotations to the number of rotations
    //required to get the total rotation S, that satisfy the minimal
    //off-diag element limit.
    cout << "N= "<< N <<" , "<< "rotations= " << number_of_rotations << endl; 
    }
}




//-----------Task 6(end)-------------
