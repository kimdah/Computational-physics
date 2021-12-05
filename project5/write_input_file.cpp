#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>


using namespace std;

void parameters(double h, double deltat, double T, double xc, double sc, int px,
                double yc, double sy, int py, double v0);

int width = 10;
ofstream ofile;

int main(int argc, char const *argv[]) {

  //ℎ, Δ𝑡, 𝑇, 𝑥𝑐, 𝜎𝑥, 𝑝𝑥, 𝑦𝑐, 𝜎𝑦, 𝑝𝑦 and 𝑣0.
  ofile.open("input.txt", ios::out | ios::trunc); // making sure file is overwritten
  ofile << setw(width) << "h" << setw(width) << "deltat" << setw(width) << "T"
        << setw(width) << "x_c" << setw(width) << "sigma_x" << setw(width) << "p_x"
        << setw(width) << "y_c" << setw(width) << "sigma_y" << setw(width) << "p_y"
        << setw(width) << "v_0" << endl;

  //(double h, double deltat, double T, double xc, double sc, int px, double yc, double sy, int py, double v0)
  parameters(0.005, 2.5e-5, 0.008, 0.25, 0.05, 200, 0.5, 0.05, 0, 0); // Task 7.1 w/o double slit
  parameters(0.005, 2.5e-5, 0.008, 0.25, 0.05, 200, 0.5, 0.10, 0, 1e10); // Task 7.3 w/double slit
  parameters(0.005, 2.5e-5, 0.002, 0.25, 0.05, 200, 0.5, 0.20, 0, 1e10); // Task 8

  ofile.close();

  return 0;
}

void parameters(double h, double deltat, double T, double xc, double sc, int px,
                double yc, double sy, int py, double v0){
  ofile << setw(width) << h << setw(width) << deltat << setw(width) << T
        << setw(width) << xc << setw(width) << sc << setw(width) << px
        << setw(width) << yc << setw(width) << sy << setw(width) << py
        << setw(width) << v0
        << endl;
}
