#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdexcept>
#include "../include/io.h"
#include "../include/lfp.h"
using namespace std;

// Arguments:
// argv[1] = path of the raster data file;
// argv[2] = index of target neuron;
// argv[3] = time range of spikes;
// argv[4] = filename of output raster file;
int main(int argc, const char* argv[]) {
  if (argc != 5) {
    throw runtime_error("wrong number of args");
  } else {
    // load data inputing arguments;
    string path = argv[1];
    int index = atoi(argv[2]);
    string range_str = argv[3];
    string::size_type pos = range_str.find_first_of(',', 0);
    double range[2];
    *range = atof(range_str.substr(0, pos).c_str());
    range_str.erase(0, pos + 1);
    *(range + 1) = atof(range_str.c_str());
    range_str = "";
    // Load target spike train;
    vector<double> spikes;
    Read1D(path, index, 0, spikes);
    // Output spike train;
    OutSpikeTrain(argv[4], spikes, range);
  }
  return 0;
}
