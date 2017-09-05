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
// 	argv[2] = path for first neuron file;
//	argv[3] = path for second neuron file
int main(int argc, const char* argv[]) {
	if (argc != 4) throw runtime_error("wrong number of args");
	// Preparing input args;
	vector<double> x;
	Read1D(argv[1], x, 0, 0);
	Print1D(argv[2], x, "app", 0);
	Read1D(argv[1], x, 1, 0);
	Print1D(argv[2], x, "app", 0);
	return 0;
}
