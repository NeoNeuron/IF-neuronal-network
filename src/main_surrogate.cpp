#include <iostream>
#include "../include/surrogate.h"
#include "../include/io.h"
using namespace std;

// compact function to calculate mutual information between multi-type signal
//	arguments:
//	argv[1] = path for inputing series;
//  argv[2] = path for outputing series;
int main(int argc, const char* argv[]) {
  if (argc != 3) throw runtime_error("wrong number of args");
  vector<double> data;
  Read1D(argv[1], data, 0, 1);
  Shuffle(data);
  Print1D(argv[2], data, "trunc", 1);
  return 0;
}
