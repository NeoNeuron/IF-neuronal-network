#include<iostream>
#include<fstream>
#include<cstdlib>
#include<iomanip>
#include<cmath>
#include"tdlc.h"
#include"../tdmi/mi.h"
using namespace std;

const double PI=3.1415926;

int main(int argc, const char* argv[]) {
  // Input args:
  //  argv[1] = timing step, with unit milliseconds;
  //	argv[2] = symetric bond of time-delay range;
  //  argv[3] = output filename;
  if (argc != 4) throw runtime_error("wrong number of args");
  // prepare data loading system;
  int tot_length = pow(2, 15);
  double dt = pow(2, -3);
  vector<double> first(tot_length), second(tot_length);
  for (int i = 0; i < tot_length; i ++) {
    first[i] = sin(i * dt);
    second[i] = sin(i * dt + PI);
  }
  // Main calculation:
  double timing_step = pow(2, -1 * atoi(argv[1]));
  int negative_time_delay = atoi(argv[2]);
  int positive_time_delay = atoi(argv[2]);
  vector<double> tdlc;
  // manipulate data;
  int n = round(timing_step / dt);
  int first_length = tot_length / n;
  vector<double> mean_first(first_length);
  double addition = 0;
  for (int i = 0; i < first_length; i ++) {
    for (int j = 0; j < n; j ++) {
      addition += first[i * n + j];
    }
    addition /= n;
    mean_first[i] = addition;
    addition = 0;
  }
  vector<double> mean_second(first_length, 0);
  addition = 0;
  for (int i = 0; i < first_length; i ++) {
    for (int j = 0; j < n; j ++) {
      addition += second[i * n + j];
    }
    addition /= n;
    mean_second[i] = addition;
    addition = 0;
  }
  // cout << mean_first.size() << '\t' << mean_second.size() << endl;
  TDLC(mean_first, mean_second, negative_time_delay, positive_time_delay, tdlc);

  // output results:
  ofstream ofile;
  ofile.open(argv[3]);
  for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i ++) {
    ofile << (double) (i - negative_time_delay) * timing_step << '\t' << setprecision(6) << (double)tdlc[i];
    if (i != (negative_time_delay + positive_time_delay)) ofile << endl;
  }
  ofile.close();
  return 0;
}
