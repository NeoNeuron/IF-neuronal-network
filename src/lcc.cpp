#include<cmath>
#include<iostream>
#include<fstream>
#include<numeric>
#include<cstdlib>
#include"../include/lcc.h"

using namespace std;

double Mean(vector<int>& x) {
  double sum = accumulate(x.begin(), x.end(), 0.0);
  return sum * 1.0 / x.size();
}

double Mean(vector<double>& x) {
  double sum = accumulate(x.begin(), x.end(), 0.0);
  return sum / x.size();
}

double Std(vector<int>& x) {
  double x_mean = Mean(x);
  vector<double> x2(x.size());
  for (int i = 0; i < x.size(); i ++) x2[i] = x[i] * x[i];
  return Mean(x2) - x_mean * x_mean;
}

double Std(vector<double>& x) {
  double x_mean = Mean(x);
  vector<double> x2(x.size());
  for (int i = 0; i < x.size(); i ++) x2[i] = x[i] * x[i];
  return Mean(x2) - x_mean * x_mean;
}

double LC(vector<int>& raster, vector<double>& lfp) {
  double raster_mean = Mean(raster);
  double lfp_mean = Mean(lfp);
  vector<double> products(raster.size());
  for (int i = 0; i < raster.size(); i ++) {
    products[i] = raster[i] * lfp[i];
  }
  return (Mean(products) - raster_mean * lfp_mean) / (Std(raster) * Std(lfp));
}

double LC(vector<double>& first, vector<double>& second) {
  double first_mean = Mean(first);
  double second_mean = Mean(second);
  vector<double> products(first.size());
  for (int i = 0; i < first.size(); i ++) {
    products[i] = first[i] * second[i];
  }
  return (Mean(products) - first_mean * second_mean) / (Std(first) * Std(second));
}

void TDLC(vector<int>& raster, vector<double>& lfp, int negative_time_delay, int positive_time_delay, vector<double>& tdlc) {
  tdlc.resize(positive_time_delay + negative_time_delay + 1);
  vector<int> raster_copy = raster;
  vector<double> lfp_copy = lfp;
  for (int i = 0; i < negative_time_delay; i++) {
    raster_copy.erase(raster_copy.begin());
    lfp_copy.erase(lfp_copy.end() - 1);
    tdlc[negative_time_delay - i - 1] = LC(raster_copy, lfp_copy);
  }
  tdlc[negative_time_delay] = LC(raster, lfp);
  raster_copy = raster;
  lfp_copy = lfp;
  for (int i = negative_time_delay + 1; i < negative_time_delay + positive_time_delay + 1; i++) {
    raster_copy.erase(raster_copy.end() - 1);
    lfp_copy.erase(lfp_copy.begin());
    tdlc[i] = LC(raster_copy, lfp_copy);
  }
}

void TDLC(vector<double>& first, vector<double>& second, int negative_time_delay, int positive_time_delay, vector<double>& tdlc) {
  tdlc.resize(positive_time_delay + negative_time_delay + 1);
  vector<double> first_copy = first;
  vector<double> second_copy = second;
  for (int i = 0; i < negative_time_delay; i++) {
    first_copy.erase(first_copy.begin());
    second_copy.erase(second_copy.end() - 1);
    tdlc[negative_time_delay - i - 1] = LC(first_copy, second_copy);
  }
  tdlc[negative_time_delay] = LC(first, second);
  first_copy = first;
  second_copy = second;
  for (int i = negative_time_delay + 1; i < negative_time_delay + positive_time_delay + 1; i++) {
    first_copy.erase(first_copy.end() - 1);
    second_copy.erase(second_copy.begin());
    tdlc[i] = LC(first_copy, second_copy);
  }
}
