#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include "../include/io.h"
using namespace std;


void OutPotential(string filename, vector<double>& v, double* t_range) {
  double dt = 0.03125;
  int tmin = floor(t_range[0] / dt);
  int tmax = floor(t_range[1] / dt);
  vector<double> v_copy(tmax - tmin);
  for (int i = tmin; i < tmax; i ++) {
    v_copy[i - tmin] = v[i];
  }
	string path = "./data/potential/" + filename;
	Print1D(path, "trunc", 1, v_copy);
}

// Arguments:
// argv[1] = path of the potential data file;
// argv[2] = index of target neuron;
// argv[3] = time range of spikes;
// argv[4] = filename of output potential file;
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
    vector<double> v;
    Read1D(path, index, 1, v);
    // Output spike train;
    OutPotential(argv[4], v, range);
  }
  return 0;
}
