#include "../include/io.h"
#include <cmath>
#include <stdexcept>
using namespace std;

int main(int argc, const char* argv[]) {
  // Input args:
  //  argv[1] = path of x;
  //  argv[2] = path of output file;
  //  argv[3] = index;
  //  argv[4] = axis;
  //  argv[5] = range, seperated by comma;
  if (argc != 6) throw runtime_error("wrong number of args");
  vector<double> data;
  Read1D(argv[1], data, atoi(argv[3]), atoi(argv[4]));
  string range = argv[5];
  int range_min = atoi(range.c_str());
  string::size_type pos = range.find_first_of(',', 0);
  int range_max = atoi(range.c_str() + pos + 1);
  double dt = 0.03125;
  data.erase(data.begin() + floor(range_min / dt), data.begin() + floor(range_max / dt));
  Print1D(argv[2], data, "trunc", atoi(argv[4]));
  return 0;
}
