#include <iostream>

using namespace std;

int main(int argc, char const *argv[]) {


  // ------ FROM EXAMPLE CODE-------
  // Check number of command line arguments
  assert(argc == 4);

  // Read command line arguments
  const int n_A = atoi(argv[1]);
  const int n_cycles = atoi(argv[2]);
  const string output_file_name = argv[3];
  // Prepare for file output
 const int print_prec = 10;
 ofstream ofile;
 ofile.open(output_file_name.c_str(), ofstream::trunc);
  /* code */
  return 0;
}
