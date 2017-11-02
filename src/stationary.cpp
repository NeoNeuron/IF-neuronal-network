#include <iostream>
#include "../include/stationary.h"
using namespace std;

void Means(vector<vector<double> > & data, vector<double>& means) {
  means.resize(data.size());
  for (size_t i = 0; i < data.size(); i ++) means[i] = Mean(data[i]);
}

void Stds(vector<vector<double> > & data, vector<double>& stds) {
  stds.resize(data.size());
  for (size_t i = 0; i < data.size(); i ++) stds[i] = Std(data[i]);
}

void Rule2(vector<vector<double> > & data, vector<vector<double> >& covs, size_t maxlag) {
  if (maxlag > data.size() - 1) maxlag = data.size() - 1;
  covs.resize(maxlag);
  double add_cov;
  for (size_t i = 0; i < maxlag - 1; i ++) {
    for (size_t j = 0; j < maxlag - i; j ++) {
        add_cov = Cov(data[j], data[j + i + 1]);
        covs[i].push_back(add_cov);
    }
  }
}

void AutoCov(vector<vector<double> >& data, vector<double>& ac, size_t indx, size_t len){
  ac.resize(len);
  vector<double> data_copy = data[indx];
  for (size_t i = 0; i < len; i ++) {
    ac[i] = Cov(data_copy, data[indx + i]);
  }
}

void AutoCov(vector<double>& data, vector<double>& ac, size_t num_lag){
  ac.resize(num_lag + 1);
  ac[0] = Cov(data, data);
  vector<double> data_copy_1 = data;
  vector<double> data_copy_2 = data;
  for (size_t i = 0; i < num_lag; i ++) {
    data_copy_1.erase(data_copy_1.begin());
    data_copy_2.erase(data_copy_2.end() - 1);
    ac[i + 1] = Cov(data_copy_1, data_copy_2);
  }
}
