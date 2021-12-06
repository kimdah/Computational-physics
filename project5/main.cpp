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
      mysstream >> h >> deltat >> T >> xc >> sx >> px >> yc >> sy >> py >> v0; //TODO: add slits to file
      if (line_counter == 0){
        // Task 7.1 w/o double slit
        problem7(Crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 0), 0);

      } else if (line_counter == 1){
        // Task 7.3 w/double slit
        int slits = 2;
        problem7(Crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, slits), slits); //commented out for testing

      } else if (line_counter == 2){
        // Task 8 and 9 (using same parameters)

        Crank crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 2); // problem 8
        problem8(crank);

        // problem9(crank, 2); // double-slit
        // problem9(Crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 1), 1);
        // problem9(Crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 3), 3);

      } else{cout << "Not all lines of the file are used" << endl;}
      line_counter +=1;
    }
  }
  else{cout << "Unable to open the file " << filename << endl;}
  myfile.close();

  return 0;
}

// --- General prob7 function to be called for 0 slits and 2 slits
void problem7(Crank crank, int slits) { // name it deviation_of_probability instead?

  cx_cube prob7 = crank.run_simulation();
  crank.to_file("U"); // what does this do?
  crank.output_probabilities(prob7, "datafiles/probability_sum_slits_"+to_string(slits)+".txt");

  // ----Remove later - only for testing:
  prob7.for_each( [](complex <double> val) { return real(conj(val)*val); } ); // changed from transform to for_each
  cube out = conv_to <cube>::from(prob7);
  out.save("datafiles/prob7_slits_" + to_string(slits) + ".dat");
  mat box = crank.make_potential_box();
  mat box_single_slit = crank.make_potential_single_slit();
  mat box_double_slit = crank.make_potential_double_slit();
  mat box_triple_slit = crank.make_potential_triple_slit();
  box.save("datafiles/box.dat");
  box_single_slit.save("datafiles/box_single_slit.dat");
  box_double_slit.save("datafiles/box_double_slit.dat");
  box_triple_slit.save("datafiles/box_triple_slit.dat");

  // ----------------
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

  cx_cube prob8 = crank.run_simulation();

  // New file with only U to get Re(u) and Im(u) : Not finished, and unsure...
  //cube out_u = conv_to <cube>::from(prob8);
  //out_u.save("datafiles/prob8_u.dat"); // other name?


  // probabilites p_{ij}^n = u_{ij}^{n*}* u_{ij}^{n} :
  prob8.for_each( [](complex <double> val) { return real(conj(val)*val); } ); // changed from transform to for_each
  cube out = conv_to <cube>::from(prob8);
  out.save("datafiles/prob8_probability.dat"); // maybe change name to something more understandable

}

void problem9(Crank crank, int slits){ // parameter slits for filenames
  /* Provides data for plotting of 1D probability function 𝑝(𝑦|𝑥=0.8;𝑡=0.002)
     for detecting the particle at screen x=0.8
     for single, double and triple slits

     Needs:
      * ?
  */
  // out.save("datafiles/prob9_slits" + to_string(slits) + ".dat");

}
