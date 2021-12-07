#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>
#include "./include/Crank.hpp"
#include <armadillo>
#include <omp.h>
#include <complex>
#include <cmath>

using namespace std::complex_literals; // to use imaginary number i | DEMANDS c++14!
using namespace std;
using namespace arma;

void problem7(Crank crank, int slits);
void problem8(Crank crank);
void problem9(Crank crank, int slits);
void workwork();

template <typename T> string to_string_with_precision(const T a_value, const int n = 2) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

int main(int argc, char const *argv[]) {
  workwork();

/*   fstream myfile;
  string filename = "input.txt";
  myfile.open(filename);
  if (myfile.is_open()){
    // This checks that the file was opened OK
    string line;
    std::getline(myfile, line); // skip the first line

    double h, deltat, T, xc, sx, px, yc, sy, py, v0, slits;
    int line_counter = 0;
    while (std::getline(myfile, line)) {
      std::stringstream mysstream(line);
      mysstream >> h >> deltat >> T >> xc >> sx >> px >> yc >> sy >> py >> v0 >> slits; //TODO: add slits to file
      if (line_counter == 0){
        // Task 7.1 w/o double slit
        
        problem7(Crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, slits), slits);

      } else if (line_counter == 1){
        // Task 7.3 w/double slit
        
        problem7(Crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, slits), slits); //commented out for testing

      } else if (line_counter == 2){
        // Task 8 and 9 (using same parameters)

        Crank crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 2); // problem 8
        problem8(crank);

        problem9(crank, 2); // double-slit
        problem9(Crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 1), 1);
        problem9(Crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, 3), 3);

      } else{cout << "Not all lines of the file are used" << endl;}
      line_counter +=1;
    }
  }
  else{cout << "Unable to open the file " << filename << endl;}
  myfile.close(); */

  return 0;
}

void workwork() {
  fstream myfile;
  string filename = "input.txt";
  myfile.open(filename);
  if (myfile.is_open()){
    // This checks that the file was opened OK
    string line;
    std::getline(myfile, line); // skip the first line
    const size_t input_vals = 15;
    std::vector<std::array<double, input_vals>> vars;
    
    int line_counter = 0;
    while (std::getline(myfile, line)) {
      double prob, h, deltat, T, xc, sx, px, yc, sy, py, v0, psum, reim, last_slice, slits;
      std::stringstream mysstream(line);
      mysstream >> prob >> h >> deltat >> T >> xc >> sx >> px >> yc >> sy >> py >> v0 >> slits >> psum >> reim >> last_slice;
      vars.push_back({{prob, h, deltat, T, xc, sx, px, yc, sy, py, v0, slits, psum, reim, last_slice}});
      line_counter +=1;
    }

    cout << std::setw( 8 ) << "Problem" << std::setw( 8 ) << "h" << std::setw( 8 ) << "deltat" << std::setw( 8 ) << "T" << std::setw( 8 ) << "x_c" << std::setw( 8 ) 
    << "sigma_x" << std::setw( 8 ) << "p_x" << std::setw( 8 ) << "y_c" << std::setw( 8 ) << "sigma_y" << std::setw( 8 ) << "p_y" << std::setw( 8 ) << "v_0" 
    << std::setw( 8 ) << "slits" << std::setw( 8 ) << "psum" << std::setw( 8 ) << "ReIm" << std::setw( 8 ) << "Last_slice" << endl; 

    for (int i = 0; i<vars.size(); i++) {
      for (int j = 0; j < (int)input_vals; j++) {
        cout << std::setw( 8 ) << vars[i][j] ;
      }
      cout << endl;
    }

    cout << "Number of simulations to run is " << line_counter << ". Running in parallel." << endl;

    #pragma omp parallel for
    for (int i = 0; i <line_counter; i++) {
      int thread_id = omp_get_thread_num();
      double prob, h, deltat, T, xc, sx, px, yc, sy, py, v0, psum, reim, last_slice, slits;
      prob=vars[i][0]; h=vars[i][1]; deltat=vars[i][2]; T=vars[i][3]; xc=vars[i][4]; sx=vars[i][5]; px=vars[i][6]; yc=vars[i][7]; sy=vars[i][8];
      py=vars[i][9]; v0=vars[i][10]; slits=vars[i][11]; psum=vars[i][12]; reim=vars[i][13]; last_slice=vars[i][14]; 

      #pragma omp critical
      cout << "Running problem " << to_string_with_precision(prob) << " on thread " << thread_id << " with parameters: h=" << h << ", delta t=" << deltat << ", T=" <<  T << ", x_c=" << xc << ", s_x=" << sx << ", p_x=" << px << ", y_c=" << yc 
      << ", s_y=" << sy << ", p_y=" << py << ", v_0=" << v0 << ", slits=" << (int)slits << ", psum=" << psum << ", ReIm=" << reim << ", last slice=" << last_slice << "." << endl;

      Crank crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, slits);
      cx_cube results_cube;
      cx_mat last_slice_mat;
      if(last_slice != 1) {
        results_cube = crank.run_simulation();
      } else if (last_slice == 1) {
        last_slice_mat = crank.run_simulation(1);
      } else {}
      // For problem 7
      if (psum == 1) {
        crank.output_probabilities(results_cube, "datafiles/Problem_"+to_string_with_precision(prob)+"_output_probability_sum_slits_"+to_string((int)slits)+".txt");
      } 
      // For problem 8-2. Results cube before being converted to probabilities. Raw Re and Im.
      if (reim == 1) {
        results_cube.save("datafiles/Problem_"+to_string_with_precision(prob)+"_ReIm_outputCube_slits_" + to_string((int)slits) + ".dat");
      }

      if(last_slice != 1) {
        results_cube.for_each( [](complex <double> val) { return real(conj(val)*val); } ); // changed from transform to for_each
        cube output = conv_to <cube>::from(results_cube);
        output.save("datafiles/Problem_"+to_string_with_precision(prob)+"_outputCube_slits_" + to_string((int)slits) + ".dat");
      } else if (last_slice == 1) {
        last_slice_mat.for_each( [](complex <double> val) { return real(conj(val)*val); } );
        mat output = conv_to <mat>::from(last_slice_mat);
        output.save("datafiles/Problem_"+to_string_with_precision(prob)+"_outputMat_slits_" + to_string((int)slits) + ".dat");
      } else {}
    }
  }
  else{cout << "Unable to open the file " << filename << endl;}
  myfile.close();
}
// // --- General prob7 function to be called for 0 slits and 2 slits
// void problem7(Crank crank, int slits) { // name it deviation_of_probability instead?

//   cx_cube prob7 = crank.run_simulation();
//   crank.to_file("U"); // what does this do?
//   crank.output_probabilities(prob7, "datafiles/probability_sum_slits_"+to_string(slits)+".txt");

//   // ----Remove later - only for testing:
//   prob7.for_each( [](complex <double> val) { return real(conj(val)*val); } ); // changed from transform to for_each
//   cube out = conv_to <cube>::from(prob7);
//   out.save("datafiles/prob7_slits_" + to_string(slits) + ".dat");
//   mat box = crank.make_potential_box();
//   mat box_single_slit = crank.make_potential_single_slit();
//   mat box_double_slit = crank.make_potential_double_slit();
//   mat box_triple_slit = crank.make_potential_triple_slit();
//   box.save("datafiles/box.dat");
//   box_single_slit.save("datafiles/box_single_slit.dat");
//   box_double_slit.save("datafiles/box_double_slit.dat");
//   box_triple_slit.save("datafiles/box_triple_slit.dat");

//   // ----------------
// }

// void problem8(Crank crank){
//   /* Provides data for colormapping time evolution of 2D probability
//      function using double slit configuration.

//      Needs:
//         *  p_{ij}^n = u_{ij}^{n*}* u_{ij}^{n}
//         *  Re(u_ij)
//         *  Im(u_ij)
//      for timesteps n corresponding to t = {0, 0.0001, 0.0002}
//   */

//   cx_cube prob8 = crank.run_simulation();

//   // New file with only U to get Re(u) and Im(u) : Not finished, and unsure...
//   //cube out_u = conv_to <cube>::from(prob8);
//   //out_u.save("datafiles/prob8_u.dat"); // other name?


//   // probabilites p_{ij}^n = u_{ij}^{n*}* u_{ij}^{n} :
//   prob8.for_each( [](complex <double> val) { return real(conj(val)*val); } ); // changed from transform to for_each
//   cube out = conv_to <cube>::from(prob8);
//   out.save("datafiles/prob8_probability.dat"); // maybe change name to something more understandable

// }

// void problem9(Crank crank, int slits){ // parameter slits for filenames
//   /* Provides data for plotting of 1D probability function 𝑝(𝑦|𝑥=0.8;𝑡=0.002)
//      for detecting the particle at screen x=0.8
//      for single, double and triple slits

//      Needs:
//       * ?
//   */
//   // out.save("datafiles/prob9_slits" + to_string(slits) + ".dat");

// }
