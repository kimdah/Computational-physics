
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <chrono>
#include <armadillo>


using namespace std;
using namespace arma;
// Performs simulations based on parameter inputs
int change_index(int i, int j, int M);
void make_matrix(int n, float r,vector<float> a);



int main(int argc, char const *argv[]) {
  int M=5;
  double deltat;
  double h;
  sp_cx_mat A= sp_cx_mat((m-2)^2, (m-2)^2).fill:zeros;
  sp_cx_mat B= sp_cx_mat((m-2)^2, (m-2)^2).fill:zeros;
  sp_cx_mat V= sp_cx_mat((m-2)^2, (m-2)^2).fill:zeros;
  cx_vec a = cx_vec((m-2)^2-1);
  cx_vec b = cx_vec((m-2)^2-1);

  for(int k=0 ; k < (m-2)^2 ; k++){
    a(k) = (1 + 4*r, deltat/2 * V(k,k);
    b(k) = (1 - 4*r, -deltat/2 * V(k,k);
    A(k,k) = a(k);
    B(k,k) = b(k);
  }

//cout << change_index(3,3,5)<<endl;



}

//changes index for the u vector(column), i and j can have values from 0 to M-2
int change_index(int i, int j, int M){return ((i%(M-1))-1)+ (M-2)*(j-1);}
