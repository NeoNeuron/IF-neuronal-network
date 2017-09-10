//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-06-03
//	Description: Mutual information analysis program;
//***************
#include "../include/mi_uniform.h"
#include "../include/io.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

using namespace std;
// compact function to calculate mutual information between multi-type signal
//	arguments:
//	argv[1] = path for double series x;
//	argv[2] = path for double series y;
//  argv[3] = number of partition of pdf of x;
//  argv[4] = number of partition of pdf of y;
int main(int argc, const char* argv[]) {
	if (argc != 5) throw runtime_error("wrong number of args");
	// Preparing input args;
	vector<double> x, y;
	Read1D(argv[1], x, 0, 1);
	Read1D(argv[2], y, 0, 1);
	size_t x_bin_num = atoi(argv[3]), y_bin_num = atoi(argv[4]);
  double mi = MI(x, y, x_bin_num, y_bin_num);
  // Output data:
  ofstream ofile;
  ofile.open("./data/mi/mi.csv", ios::app);
  ofile << mi << endl;
  ofile.close();
	return 0;
}
