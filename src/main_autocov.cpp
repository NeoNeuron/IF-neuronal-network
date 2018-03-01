#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include"../include/stationary.h"
#include"../include/io.h"
using namespace std;

int main(int argc, const char* argv[]) {
  // Input args:
  //  argv[1] = path of x;
  //  argv[2] = length of autocov;
  //  argv[3] = starting point;
  //  argv[4] = number of trials;
  if (argc != 4) throw runtime_error("wrong number of args");
  // prepare data loading system;
  vector<vector<double> > x;
  Read2D(argv[1], x);
  // Main calculation:
  size_t len = atoi(argv[2]);
  size_t indx = atoi(argv[3]);
  size_t trials = atoi(argv[4]);
  vector<double> ac(len);
  vector<double> acs(len, 0);
  for (size_t i = 0; i < trials; i ++) {
    AutoCov(x, ac, indx + i, len);
    for (size_t j = 0; j < len; j ++) acs[j] += ac[j];
  }
  for (size_t k = 0; k < len; k ++) acs[k] /= trials;
  // output results:
  Print1D("./data/stationary/autocov.csv", acs, "trunc", 1);
  return 0;
}
