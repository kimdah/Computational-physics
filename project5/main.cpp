
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <armadillo>


using namespace std;
using namespace arma;
// Performs simulations based on parameter inputs
int change_index(int i, int j, int M);
void make_matrix(int n, float r,vector<float> a);



int main(int argc, char const *argv[]) {

//cout << change_index(3,3,5)<<endl;

}

//changes index for the u vector(column), i and j can have values from 0 to M-2
int change_index(int i, int j, int M){return ((i%(M-1))-1)+ (M-2)*(j-1);}
