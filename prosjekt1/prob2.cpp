#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
using namespace std;

int main(int argc, const char * argv[]) {

    int i;
    int npoints = 10;
    double u[npoints];
    double x[npoints];
    double x_min = 0.0;
    double x_max = 1.0;
    double h = (x_max - x_min) / npoints;

    int width = 12;
    int prec = 4;

    //opening file
    ofstream ofile;
    ofile.open ("data.txt");

    //setting up the x-array and the solutions to the function u, and printing it to file
    for (i=0 ; i <= npoints ; i++){

        x[i] = h*i;
        u[i] = 1 - (1 - exp(-10)) * x[i] - exp(-10 * x[i]);
        ofile << setw(width) << setprecision(prec) << scientific << x[i] << setw(width) << setprecision(prec) << scientific << u[i] << endl;

    }
    //close file
    ofile.close();

    return 0;
}
