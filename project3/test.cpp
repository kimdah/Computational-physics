#include <armadillo>
#include <iostream>


int main(int argc, char const *argv[]) {
  /* code */
  arma::vec a = arma::vec(2).fill(2);
  arma::vec b = arma::vec(2).fill(1);


  //std::cout << a+b;
  std::cout << "norm a: " << arma::norm(a);
  std::cout << "norm b: " << arma::norm(b);
  std::cout << "norm a-b: " << arma::norm(a-b);
  //std::cout << "\n abs: " << abs(a-b);



  return 0;
}
