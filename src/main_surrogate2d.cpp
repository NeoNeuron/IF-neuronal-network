#include <iostream>
#include "../include/surrogate.h"
#include "../include/io.h"
using namespace std;

// compact function to calculate mutual information between multi-type signal
//	arguments:
//	argv[1] = path for inputing series;
//  argv[2] = path for outputing series;
//  argv[3] = axis, 0 for horizantal, 1 for vertical;
int main(int argc, const char* argv[]) {
  if (argc != 4) throw runtime_error("wrong number of args");
  vector<vector<double> > data;
  Read2D(argv[1], data);
  int axis = atoi(argv[3]);
  Shuffle2D1D(data, axis);
  Print2D(argv[2], data, "trunc");
  return 0;
}
