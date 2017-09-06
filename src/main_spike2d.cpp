#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include "../include/io.h"
#include "../include/spike.h"
using namespace std;

// Arguments:
// argv[1] = path of the rasters data file;
// argv[2] = time range of spikes; if None, then consider the whole range of spikes;
// argv[3] = binning size of binary time series of spike train;
// argv[4] = filename of output raster file;
int main(int argc, const char* argv[]) {
  if (argc != 5) throw runtime_error("wrong number of args");
  // load data inputing arguments;
  vector<vector<double> > spikes;
  Read2D(argv[1], spikes);
  string range_str = argv[2];
  double range[2];
  if (range_str == "None") {
    vector<double>::iterator p_max, p_max_current;
    p_max = max_element(spikes.begin()->begin(), spikes.end()->end());
    for (vector<vector<double> >::iterator it = spikes.begin(); it != spikes.end(); it ++) {
      p_max_current = max_element(it->begin(), it->end());
      if (*p_max_current > *p_max) p_max = p_max_current;
    }
    range[1] = ceil(*p_max);
    range[0] = 0;
  } else {
    string::size_type pos = range_str.find_first_of(',', 0);
    *range = atof(range_str.c_str());
    *(range + 1) = atof(range_str.c_str() + pos + 1);
    // Truncate the spiking series;
    for (vector<vector<double> >::iterator it = spikes.begin(); it != spikes.end(); it ++) Truncate(*it, range);
  }
  // Convert double to binary;
  vector<vector<bool> > binary_spikes(spikes.size());
  double tmax = range[1] - range[0];
  double dt = atof(argv[3]);
  for (size_t i = 0; i < binary_spikes.size(); i ++) {
    Spike2Bool(spikes[i], binary_spikes[i], tmax, dt);
  }
  // Output spike train;
  Print2D(argv[4], binary_spikes, "trunc");
  return 0;
}
