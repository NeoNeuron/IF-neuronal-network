#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include "../include/sta.h"
#include "../include/io.h"
using namespace std;

//	arguments:
//	argv[1] = path of spike train file;
//	argv[2] = path of lfp file;
//	argv[3] = tau time range, seperated by comma, with unit in milliseconds;
int main(int argc, const char* argv[]) {
  if (argc != 4) 		throw runtime_error("wrong number of args");
  // load file paths;
  string file_path_spike_train = argv[1];
  string file_path_conti_seq = argv[2];
  // load range of time delay;
  string tau_range = argv[3];
  string::size_type pos = tau_range.find_first_of(',', 0);
  double tau_min = atof(tau_range.substr(0, pos).c_str());
  tau_range.erase(0, pos + 1);
  double tau_max = atof(tau_range.c_str());
  tau_range = "";
  // input spike train and lfp;
  vector<double> spike_train, conti_seq;
  Read1D(file_path_spike_train, 0, 1, spike_train);
  Read1D(file_path_conti_seq, 0, 1, conti_seq);
  // take average of lfps;
  double sum = accumulate(conti_seq.begin(), conti_seq.end(), 0.0);
  double mean = sum / conti_seq.size();
  for (vector<double>::iterator it = conti_seq.begin(); it != conti_seq.end(); it ++) *it -= mean;

  // prepare timing systems;
  double dt_seq = 0.03125;
  int num_bin = ceil((tau_max - tau_min) / dt_seq);
  int num_bin_neg = ceil(abs(tau_min) / dt_seq);
  vector<double> sta(num_bin, 0);
  for (vector<double>::iterator it = spike_train.begin(); it != spike_train.end(); it ++) {
    int zero = floor(*it / dt_seq);
    for (int i = 0; i < num_bin; i ++) {
      if (zero - num_bin_neg + i < 0) sta[i] += 0;
      else sta[i] += conti_seq[zero - num_bin_neg + i] / spike_train.size();
    }
  }
  // output spike triggered average
  ofstream ofile;
  ofile.open("./data/sta.csv");
  ofile << "timelag,sta," << endl;
  ofile.close();
  vector<double> add(2, 0);
  vector<vector<double> > odata(num_bin, add);
  for (int i = 0; i < num_bin; i ++) {
    odata[i][0] = tau_min + i * dt_seq;
    odata[i][1] = sta[i];
  }
  Print2D("./data/sta.csv", "app", odata);
  return 0;
}
