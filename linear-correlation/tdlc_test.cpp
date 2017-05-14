#include<iostream>
#include<fstream>
#include<cstdlib>
#include<iomanip>
#include<cmath>
#include"tdlc.h"
#include"../tdmi/mi.h"
using namespace std;

const double PI=3.1415926;

int main() {
  // prepare data loading system;
  vector<double> first, second;
  int tot_length = pow(2, 15);
  double dt = pow(2, -3);
  vecotr<double> x(tot_length), y(tot_length);
  for (int i = 0; i < tot_length; i ++) {
    x[i] = Sin(i * dt);
    y[i] = Sin(i * dt + PI);
  }
  // Main calculation:
  double timing_step = pow(2, -3);
  int negative_time_delay = 50;
  int positive_time_delay = 50;
  vector<double> tdlc;
  // manipulate data;
  int n = round(timing_step / dt);
  int first_length = tot_length / n;
  vector<double> mean_first(first_length);
  double addition = 0;
  for (int i = 0; i < first_length; i ++) {
    for (int j = 0; j < n; j ++) {
      addition += first[i + j];
    }
    addition /= n;
    mean_first[i] = addition;
    addition = 0;
  }
  vector<int> mean_second(first_length, 0);
  addition = 0;
  for (int i = 0; i < second.size(); i ++) {
    for (int j = 0; j < n; j ++) {
      addition += second[i + j];
    }
    addition /= n;
    mean_second[i] = addition;
    addition = 0;
  }
  // cout << mean_first.size() << '\t' << mean_second.size() << endl;
  TDLC(mean_first, mean_second, negative_time_delay, positive_time_delay, tdlc);

  // output results:
  ofstream ofile;
  ofile.open(ofilename.c_str());
  for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i ++) {
    ofile << (double) (i - negative_time_delay) * timing_step << '\t' << setprecision(6) << (double)tdlc[i];
    if (i != (negative_time_delay + positive_time_delay)) ofile << endl;
  }
  ofile.close();
  return 0;
}
