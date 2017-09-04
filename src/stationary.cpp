#include <iostream>
#include "../include/stationary.h"
#include "../include/io.h"
using namespace std;

void Rule1(vector<vector<double> > & data, vector<double>& means) {
  means.resize(data.size());
  for (size_t i = 0; i < data.size(); i ++) {
    means[i] = Mean(data[i]);
  }
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
