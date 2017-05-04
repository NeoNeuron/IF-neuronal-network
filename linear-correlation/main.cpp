#include<iostream>
#include<vector>
#include<fstream>
#include<cstdlib>
#include<iomanip>
#include<cmath>
#include"tdlc.h"

using namespace std;

int main() {
  // prepare data loading system;
  string loading_dir = "lfp/file-txt/";
  string outputing_dir = "linear-correlation/figures/";
  string ifilename_lfp = loading_dir + "lfp.txt";
  string ifilename_raster = loading_dir + "raster.txt";
  string ofilename = outputing_dir + "linecorr.txt";
  vector<double> lfp, raster;
  ReadData(ifilename_lfp, lfp);
  ReadData(ifilename_raster, raster);
  double sampling_dt = 0.03125;
  // Main calculation:
  int negative_time_delay = 60;
  int positive_time_delay = 60;
  double timing_step = 0.25;
  vector<double> tdlc;
  // manipulate data;
  int n = round(timing_step / sampling_dt);
  int lfp_length = lfp.size() / n;
  vector<double> mean_lfp(lfp_length);
  double addition = 0;
  for (int i = 0; i < lfp_length; i ++) {
    for (int j = 0; j < n; j ++) {
      addition += lfp[i + j];
    }
    addition /= n;
    mean_lfp[i] = addition;
    addition = 0;
  }
  vector<int> mean_raster(lfp_length, 0);
  int ind;
  for (int i = 0; i < lfp_length; i ++ ) {
    ind = floor(raster[i] / timing_step);
    if (ind == lfp_length) ind -= 1;
    mean_raster[ind] = 1;
  }
  TDLC(mean_raster, mean_lfp, negative_time_delay, positive_time_delay, tdlc);

  // output results:
  ofstream ofile;
  ofile.open(ofilename.c_str());
  for (vector<double>::iterator it = tdlc.begin(); it != tdlc.end(); it ++) {
    ofile << setprecision(6) << (double)*it;
    if (it != (tdlc.end() - 1) ) ofile << endl;
  }
  ofile.close();
  return 0;
}
