#include "../include/Crank.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>

#include <armadillo>

#include <complex>
#include <cmath>

using namespace std::complex_literals; // to use imaginary number i | DEMANDS c++14!
using namespace std;
using namespace arma;

Crank::Crank(double h, double deltat, double T, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double v_0, int slits=0) {
  int M = 1/h+1; //To avvoid using M as a paramater
  t_steps_ = round(T/deltat)+1;
  M_ = M;
  h_ = h;
  deltat_ = deltat;
  poutput_ = true;
  U_empty = cx_mat(M_, M_).fill(0); // Makes a blank canvas to be reused by the col_to_mat function
  r_ = 1i*deltat/(2*pow(h,2)); //definition of r

  if(slits==0){
    V_ = make_potential_box(v_0);
  }
  if(slits==1){
    V_ = make_potential_single_slit(v_0);
  }
  if(slits==2){
    V_ = make_potential_double_slit(v_0);
  }
  if(slits==3){
    V_ = make_potential_triple_slit(v_0); // TODO: fix indexing in method
  }


  //makes matrices A and B
  make_matrices(M_, h, deltat, V_, r_); // random variables!! change

  // Commented out to test errors:
  //U_ = make_insert_wavepacket(M_, h, 0.25, 0.5, 0.05, 0.05, 200.0, 0.0); // (int M, double h, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y)
  U_ = make_insert_wavepacket(M_, h, x_c, y_c, sigma_x, sigma_y, p_x, p_y); //Tryin to manually correct syntax error where dims are swapped.

}

cx_cube Crank::run_simulation() {
  cx_vec u = construct_u_vec(U_, true); // calculates the initial u column vector
  cx_cube results = cx_cube(M_, M_, t_steps_);
  results.slice(0) = U_; // Add initial state to results cube
  cx_vec u_next;
  for (int i = 1; i<t_steps_; i++) {
    u_next = time_step(u);
    results.slice(i) = col_to_mat(u_next);
    u = u_next;
  }
  U_ = results.slice(t_steps_-1);
  return results;
}

//Problem 7
void Crank::output_probabilities(cx_cube R, string filename) {
  vec probability_sums = vec(t_steps_);
  ofstream ofile;
  ofile.open(filename);
  int width = 16;
  int prec = 4;
  ofile << setw(width) << setprecision(prec) << "Time";
  ofile << setw(width) << setprecision(prec) << "P_tot";
  ofile << endl;
  for (int i = 0; i < t_steps_; i++) {
    ofile << setw(width) << setprecision(prec) << i*deltat_;
    ofile << setw(width) << setprecision(prec) << sum_probabilies(R.slice(i));
    ofile << endl;
  }
  ofile.close();
}
// Problem 2-1
// Translates matrix (i,j) index to column (k) index which have values from 1 to M-1
int Crank::get_k_index(int i, int j, int M){
  return ((j - 1) * (M - 2)) + (i - 1);
}

cx_mat Crank::col_to_mat(cx_vec u) {
  U_empty.submat(1, 1, M_-2, M_-2) = reshape(u, M_-2, M_-2);
  return U_empty;
}

// Problem x
// Initialise potential V for the box (time-independent) - make complex if necessary
mat Crank::make_potential_box(double v0){
  mat V = mat(M_,M_); //
  V.col(0) = vec(M_).fill(v0);
  V.col(M_-1) = vec(M_).fill(v0);
  V.row(0) = rowvec(M_).fill(v0);
  V.row(M_-1) = rowvec(M_).fill(v0);
  //submat(first_row, first_col, last_row, last_col)
  V.submat(1, 1, M_-2, M_-2) = mat(M_-2,M_-2).fill(0); // filling inner matrix
  return V;
}

// Creates the potential for the double slit and box
mat Crank::make_potential_double_slit(double v0){
  mat V = make_potential_box(v0); //Creates a box of size M_* M_
  int center_index = (M_)*0.5; //200 * 0.5 = 100
  int x_thickness = 0.02/h_;// indices i in x direction: (0.02/0.005) + 1 = 5
  int x_start = center_index - x_thickness/2;
  int x_end = center_index + x_thickness/2;
  int aperture = (0.05/h_) + 1 ;// (0.05/0.005) + 1 = 11
  int center_wall_length = (0.05/h_);// (0.05/0.005) = 10. From j=95 to j=105
  int start = center_index - (center_wall_length)/2 - aperture; // j=84 for lower aperture and j=106 for upper aperture
  int end = start + aperture; // j=94 for lower aperture and j=116 for upper aperture

  for (int i = x_start; i < x_end+1; i++) {
    V.row(i).fill(v0);
    for(int j = start; j < end+1; j++) {
      V(i,j) = 0;
      V(i,j+center_wall_length+1+aperture) = 0;
    }
  }
  return V;
}

mat Crank::make_potential_single_slit(double v0){
  mat V = make_potential_box(v0); //Creates a box of size M_* M_
  int center_index = (M_)*0.5; //200 * 0.5 = 100
  int x_thickness = 0.02/h_;// indices i in x direction: (0.02/0.005) + 1 = 5
  int x_start = center_index - x_thickness/2;
  int x_end = center_index + x_thickness/2;
  int aperture = (0.05/h_) ;// (0.05/0.005) + 1 = 11
  int start = center_index - aperture/2;
  int end = start + aperture + 1;

  for (int i = x_start; i<x_end+1; i++) {
    V.row(i).fill(v0);
    for(int j = start; j < end+1; j++) {
      V(i,j) = 0;
    }
  }
  return V;
}

mat Crank::make_potential_triple_slit(double v0){
  mat V = make_potential_box(v0); //Creates a box of size M_* M_
  int center_index = (M_)*0.5; //200 * 0.5 = 100
  int x_thickness = 0.02/h_;// indices i in x direction: (0.02/0.005) + 1 = 5
  int x_start = center_index - x_thickness/2;
  int x_end = center_index + x_thickness/2;
  int aperture = (0.05/h_) +1 ;// (0.05/0.005) + 1 = 11
  int wall_length = (0.05/h_) + 1;// (0.05/0.005) + 1 = 11
  int start = center_index - wall_length - aperture - ((aperture-1)/2);
  int end = start + aperture;
  int unit_separation = wall_length + aperture;

  for (int i = x_start; i < x_end+1; i++) {
    V.row(i).fill(v0);
    for(int j = start; j < end+1; j++) {
      V(i,j) = 0;
      V(i,j + unit_separation + 1) = 0;
      V(i,j + (unit_separation * 2) + 1) = 0;
    }
  }
  return V;
}


cx_mat Crank::make_insert_wavepacket(int M, double h, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y){

  cx_mat U = cx_mat(M, M).fill(0); //Creates the matrix U. Is sparse the best choice here?

  // Find index-dimensions of wavepacket
  int x_center = (M_-1)*x_c; // 200*0.25 = index 50
  int y_center = (M_-1)*y_c; // 200*0.5 = index 100
  int x_start = round((M_-1)*x_c - ((sigma_x/h)/2)); // sigma_x is here 0.05. h is 0.005. Thus 10 h in sigma_x --> 11 points
  int x_end = x_start + (sigma_x/h) + 1;
  int y_start = round((M_-1)*y_c - ((sigma_y/h)/2));
  int y_end = y_start + (sigma_y/h) + 1;

  // If wavepacket starts overlapping or outside a boundary, try to relocate it inside the boundary
  if (x_start < 1) {
    cout << "Wavepacket too far to the left, attempting to move to the right." << endl;
    int move_right = -x_start + 1;
    x_start = x_start + move_right;
    x_end = x_end + move_right;
  }
  if (y_start < 0) {
    cout << "Wavepacket too far to the up, attempting to move down." << endl;
    int move_down = -y_start + 1;
    y_start = y_start + move_down;
    y_end = y_end+ move_down;
  }
  if (x_end > M_-2) {
    cout << "Wavepacket too far to the right, attempting to move to the left." << endl;
    int move_left = x_end - M_ +2;
    x_end = x_end - move_left;
    x_start = x_start - move_left;
  }

  if (y_end > M_-2) {
    cout << "Wavepacket too far to the down, attempting to move up." << endl;
    int move_up = y_end - M_ - 2;
    y_end = y_end - move_up;
    y_start = y_start - move_up;

  }
  // cout << endl;
  // cout << "x_start:" << x_start << " x_end: " << x_end << endl;
  // cout << "y_start:" << y_start << " y_end: " << y_end << endl;

  //Then check that it is still within boundary, if else then it's too big
  if (x_start < 1) {cout << "Too far left " << endl;}
  if (y_start < 1) {cout << "Too far up " << endl;}
  if (x_end > M-2) {cout << "Too far right " << endl;}
  if (y_end > M-2) {cout << "Too far down " << endl;}
  if (x_start < 1 || y_start < 1 || x_end > M-2 || y_end > M-2) {cout << "That's what she said" << endl;}
  double psum = 0;

  // Inserts the wavepacket and calculates normalisation factor
  x_start = 1;
  x_end = M-2;
  y_start = 1;
  y_end = M-2;
  for(int i = x_start; i< x_end; i++){
    double x = i*h;
    for(int j = y_start; j< y_end; j++){
      double y = j*h;
      complex <double> c = exp( -( pow(x - x_c ,2) / ( 2 * pow(sigma_x,2) ) ) - ( pow(y - y_c, 2) / ( 2*pow(sigma_y, 2) ) ) + 1i * p_x * (x - x_c) + 1i * p_y * (y - y_c) );
      U(i,j) = c;
      psum += real(conj(c)*c);
    }
  }
  cout << "sum of real magnitudes is " << psum << endl; // remove later
  // Normalises U to 1
  double psum2 = 0;
  for(int i = x_start; i< x_end; i++){
    for(int j = y_start; j< y_end; j++){
      complex <double> c = cx_double(1/sqrt(psum)) * cx_double(U(i,j));
      U(i,j) = c;
      psum2 += real(conj(c)*c);
    }
  }
  cout << "sum of real magnitudes is normalised to " << psum2 << endl;
  return U;
}



double Crank::sum_probabilies(cx_mat U) {
  double psum = 0;
  for(int i = 0; i< M_; i++){
    for(int j = 0; j< M_; j++){
      psum += real(conj(U(i,j))*U(i,j));
    }
  }
  return psum;
}



// constructs the u vector based on U matrix
cx_vec Crank::construct_u_vec(cx_mat U, bool normalise){
  int M = sqrt(U.size()); //size() gives M², assumes U to be quadratic
  cx_vec u = cx_vec(pow(M-2,2));

  for(int i = 1; i< M-1; i++){
    for(int j = 1; j< M-1; j++){
      u(get_k_index(i, j, M)) = U(i, j);
    }
  }
  /* if(normalise){ //Normalizes u, but not sure if correct
    sp_cx_mat normalisation_factor;
    normalisation_factor = u.t()*u; //this is a 1x1 matrix
    u = u/normalisation_factor(0, 0);
  } */
  return u;
}

// Task3: Returns the column vector of the next time step u(n+1)
cx_vec Crank::time_step(cx_vec u) {
  cx_vec b = B_ * u; // Solves eq 26 RHS with u(n)
  return spsolve(A_, b); // Solves eq 26 for u(n+1). Here the vector u(n) (and similarly u(n+1) is a column vector that contains the u^n_ij values for all the internal points of the xy grid at time step n
 }


// Makes specialized A and B matrices (Task 2.3)
void Crank::make_matrices(int M, double h, double deltat, mat V, complex<double> r){
  int mat_size = pow(M-2,2);
  cx_vec a = cx_vec(mat_size);
  cx_vec b = cx_vec(mat_size);

  for (int i = 1; i < M-1; i++){ // Excluding boundaries in V (infinity) - is that ok?
    for (int j = 1; j < M-1; j++){
      int k = get_k_index(i,j,M);
      a(k) = (1.0 + 4.0*r + 1.0i*(deltat/2*cx_double(V(i,j))));
      b(k) = (1.0 - 4.0*r - 1.0i*(deltat/2*cx_double(V(i,j))));
    }
  }
  A_ = make_matrix(-r, a);
  B_ = make_matrix(r,  b);
}

// Makes a matrix on the general form of A and B
sp_cx_mat Crank::make_matrix(complex<double> r, cx_vec d){
  int S = d.size(); // (M-2)^2
  int s = sqrt(S); // (M-2)
  sp_cx_mat M = sp_cx_mat(S, S);

  // Making diag
  for (int i = 0; i < S; i+=s){ // last index = S-s
    sp_cx_mat D(s, s);
    D.diag() = d.subvec(i, i+s-1);// diagonal
    D.diag(-1) = cx_vec(s-1).fill(r); // subdiagonal
    D.diag(1) = cx_vec(s-1).fill(r); // superdiagonal
    //submat(first_row, first_col, last_row, last_col)
    M.submat(i, i, i+s-1, i+s-1) = D; //ex: (0,0,2,2), (3,3,5,5)
  }

  // Making non-diag, non-corners
  for (int i = s; i < S; i+=s){ // last index= S-s
    sp_cx_mat ND(s,s);
    ND.diag() = cx_vec(s).fill(r); // fill diagonal with r value
    //submat(first_row, first_col, last_row, last_col)
    M.submat(i-s, i, i-1, i+s-1) = ND; //ex: i=s=3: (0,3,2,5) i=2s=6: (3,6,5,8)
    M.submat(i, i-s, i+s-1, i-1) = ND; //ex: i=s=3: (3,0,5,2) i=2s=6: (6,3,8,5)
  }

  // Corners are 0 matrix, which they are already initialized as through sp_ (sparse)
  return M;
}

void Crank::print() {
  cout << U_.size()<<endl;
  for (int i = 0; i<M_; i++) {
    for (int j = 0; j<M_; j++) {
      cout << setw(25) << U_(i,j);
    }
  cout << endl;
  }
}


// Just a function to see what a sparse matrix looks like
int Crank::to_file(string s) {
  sp_cx_mat A;
  string filename;
  int largeness;
  int width = 24;
  int prec = 3;
  if (s == "A") {
    A = A_;
    largeness = A.size();
    filename = "datafiles/sparse_matrix_A_size_" + to_string(largeness) + ".txt";
  } else if (s == "B") {
    A = B_;
    largeness = A.size();
    filename = "datafiles/sparse_matrix_B_size_" + to_string(largeness) + ".txt";
  } else if (s == "U") {
    A = U_;
    largeness = A.size();
    filename = "datafiles/sparse_matrix_U_size_" + to_string(largeness) + ".txt";
  } else {
    cout << "Unknown matrix\n";
    return 1;
  }
  largeness = sqrt(largeness);

  //A.save("datafiles/test.dat");

  ofstream ofile;
  ofile.open(filename);
  for (int i = 0; i<largeness; i++) {
    ofile << setw(width) << setprecision(prec) << i;
  }
  ofile << endl;
  for (int i = 0; i<largeness; i++) {
    for (int j = 0; j<largeness; j++) {
      ofile << setw(width) << setprecision(prec) << scientific << A(i,j);
    }
    ofile << endl;
  }
  ofile.close();
  return 0;
}


// // Creates the potential for the double slit and box
// mat Crank::make_potential_double_slit(double v0){
//   mat V = make_potential_box(v0); //Creates a box of size M_* M_

//   double h = 1.0/(M_-1);

//   //Finds the right indeces according to the dimesions specified
//   int wall_thickness_index = floor(0.02/h)/2; //0.02
//   int wall_position_index = floor(0.5/h);      //0.5
//   int slit_seperation_index = floor(0.05/h)/2;//0.05
//   int slit_epeture_index = floor(0.05/h);      //0.05

//   //Sets up how the wallshould look
//   vec wall_config = vec(M_).fill(v0);
//   for(int i=0; i<slit_epeture_index; i++){
//     wall_config(floor(0.5/h)+slit_seperation_index+i) = 0;
//     wall_config(floor(0.5/h)-slit_seperation_index-i) = 0;
//   }

//   // the wall appears to be made along the column direction. It should be row direction
//   // Builds the wall
//   // for(int i =0; i < wall_thickness_index; i++){
//   //   V.col(wall_position_index + i) = wall_config;
//   //   V.col(wall_position_index - i) = wall_config;
//   // }

//    for(int i =0; i < wall_thickness_index; i++){
//     V.row(wall_position_index + i) = wall_config;
//     V.row(wall_position_index - i) = wall_config;

//   }
//   return V;
// }

// mat Crank::make_potential_single_slit(double v0){
//   mat V = make_potential_box(v0); //Creates the box potential

//   double h = 1.0/(M_-1);

//   //Finds the righ indeces according to the dimesions specified
//   int wall_thickness_index = floor(0.02/h)/2; //0.02
//   int wall_position_index = floor(0.5/h);      //0.5
//   int slit_apeture_index = floor(0.05/h);      //0.05

//   //Sets up how the wallshould look
//   vec wall_config = vec(M_).fill(v0);
//   for(int i=0; i < slit_apeture_index; i++){
//     wall_config(floor(0.5/h)-floor(slit_apeture_index/2)+i) = 0; //floor(slit_apeture_index/2) -->centering the slit
//   }
//   //Builds the wall
//   for(int i =0; i < wall_thickness_index; i++){
//     V.col(wall_position_index + i) = wall_config;
//     V.col(wall_position_index - i) = wall_config;

//   }
//   return V;
// }

// mat Crank::make_potential_triple_slit(double v0){
//   mat V = make_potential_box(v0); //Creates the box potetnial

//   double h = 1.0/(M_-1);

//   //Finds the righ indeces according to the dimesions specified
//   int wall_thickness_index = floor(0.02/h)/2; //0.02
//   int wall_position_index = floor(0.5/h);      //0.5
//   int slit_seperation_index = floor(0.05/h)/2;//0.05
//   int slit_apeture_index = floor(0.05/h);      //0.05

//   //Sets up how the wallshould look
//   vec wall_config = vec(M_).fill(v0);
//   for(int i=0; i<slit_apeture_index; i++){
//     wall_config(floor(0.5/h)+slit_seperation_index+slit_apeture_index/2+i) = 0; //lower slit
//     wall_config(floor(0.5/h)-slit_seperation_index+slit_apeture_index/2-i) = 0;//upper slit
//     wall_config(floor(0.5/h)-floor(slit_apeture_index/2)+i) = 0;              //middle slit
//   }
//   //Builds the wall
//   for(int i =0; i<wall_thickness_index;i++){
//     V.col(wall_position_index + i) = wall_config;
//     V.col(wall_position_index - i) = wall_config;

//   }
//   return V;
// }