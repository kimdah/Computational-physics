#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>
#include "./include/Crank.hpp"
#include <armadillo>

#include <complex>
#include <cmath>

using namespace std::complex_literals; // to use imaginary number i | DEMANDS c++14!
using namespace std;
using namespace arma;

void problem7();
int main(int argc, char const *argv[]) {
  // Crank(double h, double deltat)
  // ----- Cecilie-----
  /*
  fstream myfile;
  string filename = "input.txt";
  myfile.open(filename);
  if (myfile.is_open()){
    // This checks that the file was opened OK
    string line;
    std::getline(myfile, line); // skip the first line

    double h, deltat, T, xc, sc, px, yc, sy, py, v0;
    int line_counter = 0;
    while (std::getline(myfile, line)) {
      std::stringstream mysstream(line);
      mysstream << h << deltat << T << xc<< sc<< px << yc << sy << py << v0;
      //Crank crank(h, deltat, T, xc, sc, px, yc, sy, py, v0); // comment out once Crank constructor is fixed

      if (line_counter == 0){
        // Task 7.1 w/o double slit
        //problem7(crank); // send in crank-object as parameter?
        //cout << line << endl; // tester!
      } else if (line_counter == 1){
        // Task 7.3 w/double slit
      } else if (line_counter == 2){
        // Task 8
      } else{
        cout << "Not all lines of the file are used" << endl;
      }

    }

  }
  else
  {
    cout << "Unable to open the file " << filename << endl;
  }
  myfile.close();
  */
  // ---- end Cecilie-----

  problem7();

  return 0;
}

void problem7() {
  //Crank crank(0.005, 2.5e-5);
  double v_0 = numeric_limits<double>::max(); //Large potential
  //Crank crank(0.005, 2.5e-5, 0.008, 0.25, 0.5, 0.05, 0.05, 200, 0.0, v_0, 1); this should be correct
  Crank crank(0.005, 2.5e-5, 0.008, 0.5, 0.25, 0.05, 0.05, 0.0, -200.0, v_0, 1);  //Manually handling the index swap


  //crank.to_file("A");
  //crank.to_file("B");

  cx_cube prob7_1 = crank.run_simulation(321); // TODO: 321 gives T = 0.008. Change input to actual time and modify code in Crank
  crank.to_file("U");

  crank.output_probabilities(prob7_1, "datafiles/probability_sum_test.txt");
  prob7_1.transform( [](complex <double> val) { return real(conj(val)*val); } );
  cube out = conv_to< cube >::from(prob7_1);
  out.save("datafiles/prob7_1.dat");
}
