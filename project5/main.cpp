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

void problem7(Crank crank, int slits);
void problem8(Crank crank);
void problem9(Crank crank, int slits);

int main(int argc, char const *argv[]) {
  // ----- Cecilie-----

  fstream myfile;
  string filename = "input.txt";
  myfile.open(filename);
  if (myfile.is_open()){
    // This checks that the file was opened OK
    string line;
    std::getline(myfile, line); // skip the first line

    double h, deltat, T, xc, sx, px, yc, sy, py, v0;
    int line_counter = 0;
    while (std::getline(myfile, line)) {
      std::stringstream mysstream(line);
      mysstream >> h >> deltat >> T >> xc >> sx >> px >> yc >> sy >> py >> v0;
      if (line_counter == 0){
        // Task 7.1 w/o double slit

        //Crank crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 0);
        problem7(Crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 0), 0);
      } else if (line_counter == 1){
        // Task 7.3 w/double slit

        //Crank crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 2);
        problem7(Crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 2), 2); //commented out for testing
      } else if (line_counter == 2){
        // Task 8 and 9 (using same parameters)

        // Crank crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 2); // problem 8

        // problem9(Crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 1), 1);
        // problem9(Crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 2), 2);
        // problem9(Crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 3), 3);

      } else{cout << "Not all lines of the file are used" << endl;}
      line_counter +=1;
    }
  }
  else{cout << "Unable to open the file " << filename << endl;}
  myfile.close();

  // ---- end Cecilie-----

  return 0;
}

// --- General prob7 function to be called for 0 slits and 2 slits
void problem7(Crank crank, int slits) { // name it deviation_of_probability instead?
  // slits as input for filename
  double v_0 = numeric_limits<double>::max(); //Large potential - where is this used now?
  //Crank crank(0.005, 2.5e-5, 0.008, 0.25, 0.5, 0.05, 0.05, 200, 0.0, 1e10, 2); //this should be correct - shouldn't s_y be 0.1? (see 7.3)
  //Crank crank(0.005, 2.5e-5, 0.008, 0.5, 0.25, 0.05, 0.05, 0.0, 200.0, 1e10, 2);  //Manually handling the index swap

  //crank.to_file("A");
  //crank.to_file("B");

  cx_cube prob7 = crank.run_simulation();
  crank.to_file("U");

  crank.output_probabilities(prob7, "datafiles/probability_sum_slits"+to_string(slits)+".txt");
  prob7.for_each( [](complex <double> val) { return real(conj(val)*val); } ); // changed from transform to for_each
  cube out = conv_to <cube>::from(prob7);
  //out.save("plotfiles/one.dat"); // is this a test name?

  out.save("datafiles/prob7_slits_" + to_string(slits) + ".dat");
}

void problem8(Crank crank){
  /* Provides data for colormapping time evolution of 2D probability
     function using double slit configuration.

     Needs:
        *  p_{ij}^n = u_{ij}^{n*}* u_{ij}^{n}
        *  Re(u_ij)
        *  Im(u_ij)
     for timesteps n corresponding to t = {0, 0.0001, 0.0002}
  */
  

}
void problem9(Crank crank, int slits){ // parameter slits for filenames
  /* Provides data for plotting of 1D probability function 𝑝(𝑦|𝑥=0.8;𝑡=0.002)
     for detecting the particle at screen x=0.8
     for single, double and triple slits

     Needs:
      *
  */

}
