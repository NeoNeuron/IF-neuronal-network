#include "../include/io.h"
#include <stdexcept>
using namespace std;

int main(int argc, const char* argv[]) {
  // Input args:
  //  argv[1] = path of x;
  //  argv[2] = path of output file;
  //  argv[3] = index;
  //  argv[4] = axis;
  if (argc != 5) throw runtime_error("wrong number of args");
  vector<double> data;
  Read1D(argv[1], data, atoi(argv[3]), atoi(argv[4]));
  Print1D(argv[2], data, "trunc", atoi(argv[4]));
  return 0;
}
