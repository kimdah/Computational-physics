#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <armadillo>
#include <cmath>

using namespace std;
// Here we define various functions called by the main program
// this function defines the function to integrate
double func(double x);
// Main function begins here
int main(){
    int i, n;
    long idum;
    double crude_mc, crude_mc2, x, sum_sigma, fx, variance;
    cout << "Read in the number of Monte-Carlo samples" << endl;
    cin >> n;
    
    crude_mc = sum_sigma=0. ; idum=-1 ;
    // evaluate the integral with the a crude Monte-Carlo method
    for ( i = 1; i <= n; i++){
    x=((double)rand()/(double)RAND_MAX);
    fx=func(x);
    crude_mc += fx;
    sum_sigma += fx*fx;
    }

    crude_mc = crude_mc/((double) n );
    sum_sigma = sum_sigma/((double) n );
    variance=sum_sigma-crude_mc2*crude_mc2;
    // final output
    cout << " variance= " << variance << " Integral = "
    << crude_mc << " Exact= " << M_PI << endl;
} // end of main program
  // this function defines the function to integrate
double func(double x) {
    double value;
    value = 4/(1.+x*x);
    return value;
} // end of function to evaluate