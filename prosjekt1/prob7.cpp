#include <armadillo>
#include <iostream>
#include <string>


int main(int argc, char* argv[]){

  if (argc != 2){
    std::string executable_name = argv[0];

    std::cerr << "Error: Wrong number of input arguments." << std::endl;
    std::cerr << "Usage: " << executable_name << " <number of iterations (integer)> " << std::endl;

    // Exit program with non-zero return code to indicate a problem
    return 1;
  }
  int npoints = atoi(argv[1]);
  //int npoints = 10;

  double x_min = 0.0;
  double x_max = 1.0;
  double h = (x_max - x_min) / npoints;

  arma::vec a = arma::vec(n).fill(-1.);
  arma::vec b = arma::vec(n).fill(2.);
  arma::vec c = arma::vec(n).fill(-1);

  a(0) = 0.;
  c(npoints) = 0.;

  arma::vec v = arma::vec(n);

  arma::vec g = arma::vec(n);
  g(0) = h*h*100*exp(-10*(0+h)) + 0; //?


  arma::vec btilde = arma::vec(n);
  arma::vec gtilde = arma::vec(n);



  std::cout << "x0: " << x(0)<< "\n";
  return 0;

}
