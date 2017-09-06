//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-06-03
//	Description: Mutual information analysis program;
//***************
#include "../include/io.h"
#include <iostream>
#include <stdexcept>
using namespace std;
// compact function to operate stationary test;
//	arguments:
//	argv[1] = path for sampled time series x;
//	argv[2] = path for output file;
int main(int argc, const char* argv[]) {
	if (argc != 3) throw runtime_error("wrong number of args");
	// Preparing input args;
	vector<vector<double> > x;
	Read2D(argv[1], x);
	vector<vector<double> > y;
	Transpose(x, y);
	Print2D(argv[2], y, "trunc");
	return 0;
}
