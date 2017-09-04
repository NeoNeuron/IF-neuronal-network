//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-06-03
//	Description: Mutual information analysis program;
//***************
#include "../include/stationary.h"
#include "../include/io.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

using namespace std;
// compact function to operate stationary test;
//	arguments:
//	argv[1] = path for sampled time series x;
int main(int argc, const char* argv[]) {
	if (argc != 2) throw runtime_error("wrong number of args");
	// Preparing input args;
	vector<vector<double> > x;
	Read2D(argv[1], x);
	vector<double> means;
	vector<vector<double> > covs;
	Rule1(x, means);
	Print1D("./data/stationary/stat_means.csv", means, "trunc", 1);
	Rule2(x, covs, 100);
  Print2D("./data/stationary/stat_covs.csv", covs, "trunc");
	return 0;
}
