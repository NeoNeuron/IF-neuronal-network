#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <cstdlib>
#include"../include/cc.h"

using namespace std;

template <class T> double Mean(vector<T>& x) {
  double sum = accumulate(x.begin(), x.end(), 0.0);
  return sum * 1.0 / x.size();
}

template <class T> double Std(vector<T>& x) {
  double x_mean = Mean(x);
  double x2 = inner_product(x.begin(), x.end(), x.begin(), 0.0);
  return x2 / x.size() - x_mean * x_mean;
}

double PCC(vector<double>& x, vector<double>& y) {
  double product = inner_product(x.begin(), x.end(), y.begin(), 0.0);
  return (product / x.size() - Mean(x) * Mean(y)) / (Std(x) * Std(y));
}

void Correlation(vector<double>& x, vector<vector<double> >& y, vector<double>& cc) {
  cc.resize(y.size());
  for (int i = 0; i < y.size(); i++) {
    cc[i] = PCC(x, y[i]);
  }
}
