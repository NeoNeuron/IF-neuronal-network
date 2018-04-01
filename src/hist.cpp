#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include "../include/io.h"
using namespace std;

double Max(vector<double>& data) {
	return *(max_element(data.begin(), data.end()));
}

double Min(vector<double>& data) {
	return *(min_element(data.begin(), data.end()));
}

int main(int argc, const char* argv[]) {
  if (argc != 2) throw runtime_error("wrong number of args");
  // Read the continues variable in;
  vector<double> data;
  Read1DBin(argv[1], data, 0, 0);
  size_t num_pairs = data.size();
  double dmax = Max(data), dmin = Min(data);
  printf(">> Data ranging from [%.4f, %.4f]\n", dmin, dmax);
  size_t num_bin = 20;
  double binsize = (dmax - dmin) / num_bin;
  vector<double> pdf(num_bin);
  double increament = 1.0 / num_pairs;
  size_t ind;
	for (size_t i = 0; i < num_pairs; i++) {
		ind = floor((data[i] - dmin) / binsize);
		if (ind == num_bin) ind = num_bin - 1;
		pdf[ind] += increament;
	}
  Print1D("./data/tmp/pdf.csv", pdf, "trunc", 1);
  return 0;
}
