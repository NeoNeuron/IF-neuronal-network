#include <iostream>
#include <stdexcept>
#include"../include/stationary.h"
#include"../include/io.h"
using namespace std;

int main(int argc, const char* argv[]) {
  // Input args:
  //  argv[1] = path of x;
  if (argc != 2) throw runtime_error("wrong number of args");
  // prepare data loading system;
  vector<vector<double> > x;
  Read2D(argv[1], x);
  // Main calculation:
  vector<double> ac;
  AutoCov(x, ac);
  // output results:
  Print1D("./data/stationary/autocov.csv", ac, "trunc", 1);
  return 0;
}
