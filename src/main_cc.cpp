#include <iostream>
#include <stdexcept>
#include"../include/cc.h"
#include"../include/io.h"
using namespace std;

int main(int argc, const char* argv[]) {
  // Input args:
  //  argv[1] = path of x;
  //  argv[2] = path of y;
  if (argc != 3) throw runtime_error("wrong number of args");
  // prepare data loading system;
  vector<double> x;
  vector<vector<double> > y;
  Read1D(argv[1], 0, 1, x);
  Read2D(argv[2], y);
  // Main calculation:
  vector<double> cc;
  Correlation(x, y, cc);
  // output results:
  Print1D("./data/cc/cc.csv", "trunc", 1, cc);
  return 0;
}
