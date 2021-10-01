#include <armadillo>
#include <iostream>


int main(int argc, char const *argv[]) {
  /* code */
  arma::vec a = arma::vec(2).fill(2);
  arma::vec b = arma::vec(2).fill(1);


  std::cout << a+b;



  return 0;
}
