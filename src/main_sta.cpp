#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include "../include/io.h"
using namespace std;

//	arguments:
//	argv[1] = path of spike train file;
//	argv[2] = path of lfp file;
//	argv[3] = time range, seperated by comma;
int main(int argc, const char* argv[]) {
  if (argc != 4) throw runtime_error("wrong number of args");
  // input spike train and lfp;
  vector<bool> spike_train;
  vector<double> conti_seq;
  Read1D(argv[1], spike_train, 0, 1);
  Read1D(argv[2], conti_seq, 0, 1);
  // load range of time delay;
  string tau_range = argv[3];
  string::size_type pos = tau_range.find_first_of(',', 0);
  double ntd = atof(tau_range.substr(0, pos).c_str());
  tau_range.erase(0, pos + 1);
  double ptd = atof(tau_range.c_str());
  tau_range = "";
  vector<double> sta(ntd + ptd + 1, 0);
  // calculate number of spikes;
  int number = 0;
	for (vector<bool>::iterator it = spike_train.begin(); it != spike_train.end(); it ++) {
		if (*it) number ++;
	}
  if (number != 0) {
    // take average of lfps;
    double sum = accumulate(conti_seq.begin(), conti_seq.end(), 0.0);
    double mean = sum / conti_seq.size();
    for (vector<double>::iterator it = conti_seq.begin(); it != conti_seq.end(); it ++) *it -= mean;
    double add_sta;
    // No shift;
		add_sta = Product(spike_train, conti_seq);
    sta[ntd] = add_sta / number;

		// Negative shift;
		vector<bool> bool_copy = spike_train;
		vector<double> double_copy = conti_seq;
		for (int i = 0; i < ntd; i++) {
			bool_copy.erase(bool_copy.begin());
			double_copy.erase(double_copy.end() - 1);
			add_sta = Product(bool_copy, double_copy);
      sta[ntd - i - 1] = add_sta / number;
		}

		// Positive shift;
		bool_copy = spike_train;
		double_copy = conti_seq;
		for (int i = 0; i < ptd; i++) {
			bool_copy.erase(bool_copy.end() - 1);
			double_copy.erase(double_copy.begin());
			add_sta = Product(bool_copy, double_copy);
      sta[ntd + i + 1] = add_sta / number;
		}
  }
  // output spike triggered average
  ofstream ofile;
  ofile.open("./data/sta.csv");
  ofile << "timelag,sta," << endl;
  ofile.close();
  vector<vector<double> > odata(ntd + ptd + 1, vector<double>(2, 0));
  for (int i = 0; i < ntd + ptd + 1; i ++) {
    odata[i][0] = -ntd + i;
    odata[i][1] = sta[i];
  }
  Print2D("./data/sta.csv", odata, "app");
  return 0;
}
