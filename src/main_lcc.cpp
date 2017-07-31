#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include"../include/lcc.h"
#include"../include/io.h"
using namespace std;

int main(int argc, const char* argv[]) {
  // Input args:
  //  argv[1] = path of raster;
  //  argv[2] = path of lfp;
  //  argv[3] = timing step, with unit milliseconds;
  //	argv[4] = time-delay range, seperated by comma;
  if (argc != 5) throw runtime_error("wrong number of args");
  // prepare data loading system;
  string ipath_raster = argv[1];
  string ipath_lfp = argv[2];
  string opath = "./data/lcc/lcc.csv";
  vector<double> lfp, raster;
  Read1D(ipath_lfp, 0, 1, lfp);
  Read1D(ipath_raster, 0, 1, raster);
  double sampling_dt = 0.03125;
  // Main calculation:
  double timing_step = atof(argv[3]);
  string range = argv[4];
	string::size_type pos = range.find_first_of(',', 0 );
	int negative_time_delay = atoi(range.substr(0, pos).c_str());
	range.erase(0, pos + 1);
	int positive_time_delay = atoi(range.c_str());
	range = "";
  vector<double> tdlc;
  // manipulate data;
  int n = round(timing_step / sampling_dt);
  int lfp_length = lfp.size() / n;
  vector<double> mean_lfp(lfp_length);
  double addition = 0;
  for (int i = 0; i < lfp_length; i ++) {
    for (int j = 0; j < n; j ++) {
      addition += lfp[i * n + j];
    }
    addition /= n;
    mean_lfp[i] = addition;
    addition = 0;
  }
  vector<int> mean_raster(lfp_length, 0);
  int ind;
  for (int i = 0; i < raster.size(); i ++) {
    ind = floor(raster[i] / timing_step);
    if (ind == lfp_length) ind -= 1;
    mean_raster[ind] += 1;
  }
  // cout << mean_lfp.size() << ',' << mean_raster.size() << endl;
  TDLC(mean_raster, mean_lfp, negative_time_delay, positive_time_delay, tdlc);

  // output results:
  ofstream ofile;
  ofile.open(opath.c_str());
  ofile << "timelag,lcc" << endl;
  for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i ++) {
    ofile << (double) (i - negative_time_delay) * timing_step << ',' << setprecision(6) << (double)tdlc[i];
    if (i != (negative_time_delay + positive_time_delay)) ofile << endl;
  }
  ofile.close();
  return 0;
}
