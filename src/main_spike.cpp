#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <stdexcept>
#include "../include/io.h"
#include "../include/spike.h"
using namespace std;
int myrandom(int i) {return rand()%i;};
// Arguments:
// argv[1] = path of the raster data file;
// argv[2] = path of output raster file;
// argv[3] = index of target neuron;
// argv[4] = time range of spikes; if None, then consider the whole range of spikes;
// argv[5] = binning size of binary time series of spike train;
// argv[6] = shuffle flag;
int main(int argc, const char* argv[]) {
  if (argc != 7) throw runtime_error("wrong number of args");
  // load data inputing arguments;
  vector<double> spikes;
  Read1D(argv[1], spikes, atoi(argv[3]), 0);
  string range_str = argv[4];
  double range[2];
  stringstream irange(argv[4]);
  string buffer;
  for (size_t i = 0; i < 2; i++) {
    getline(irange, buffer, ',');
    range[i] = atof(buffer.c_str());
  }
  // Truncate the spiking series;
  Truncate(spikes, range);
  // Convert double to binary;
  vector<bool> binary_spikes;
  double tmax = range[1] - range[0];
  double dt = atof(argv[5]);
  Spike2Bool(spikes, binary_spikes, tmax, dt);
  // shuffle flag;
  bool shuffle_flag;
  istringstream s_in(argv[6]);
  s_in >> boolalpha >> shuffle_flag;
  srand(unsigned(time(0)));
  if (shuffle_flag) random_shuffle(binary_spikes.begin(), binary_spikes.end(), myrandom);
  // Output spike train;
  Print1D(argv[2], binary_spikes, "trunc", 1);
  return 0;
}
