//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-06-03
//	Description: Mutual information analysis program;
//***************
#include "../include/mi.h"
#include "../include/io.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

using namespace std;
// compact function to calculate mutual information between multi-type signal
//	arguments:
//	argv[1] = path for double series x;
//	argv[2] = path for double series y;
int main(int argc, const char* argv[]) {
	if (argc != 3) throw runtime_error("wrong number of args");
	// Preparing input args;
	vector<double> x, y;
	Read1D(argv[1], 0, 1, x);
	Read1D(argv[2], 0, 1, y);

  double mi = MI(x, y);
  // Output data:
  ofstream ofile;
  ofile.open("./data/mi/mi.csv", ios::app);
  ofile << mi << endl;
  ofile.close();
	return 0;
}
