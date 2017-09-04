#include <iostream>
#include <stdexcept>
#include"../include/stationary.h"
#include"../include/io.h"
using namespace std;

int main(int argc, const char* argv[]) {
  // Input args:
  //  argv[1] = path of x;
  //  argv[2] = number of time lag;
  if (argc != 3) throw runtime_error("wrong number of args");
  // prepare data loading system;
  vector<double> x;
  Read1D(argv[1], x, 0, 1);
  size_t num_lag = atoi(argv[2]);
  // Main calculation:
  vector<double> ac;
  AC(x, ac, num_lag);
  // output results:
  Print1D("./data/stationary/ac.csv", ac, "trunc", 1);
  return 0;
}
