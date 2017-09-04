#ifndef _IFNET_STATIONARY_H_
#define _IFNET_STATIONARY_H_

#include <vector>
#include <string>
#include <numeric>
using namespace std;

template <class T> double Mean(vector<T>& x){
  return accumulate(x.begin(), x.end(), 0.0) / x.size();
}

template <class T> double Std(vector<T>& x) {
  double x2 = inner_product(x.begin(), x.end(), x.begin(), 0.0);
  return x2 / x.size() - Mean(x) * Mean(x);
}

// Covariance:
//  Return: value of covariance
template <class T> double Cov(vector<T>& x, vector<T>& y) {
  double product = inner_product(x.begin(), x.end(), y.begin(), 0.0);
  return (product / x.size() - Mean(x) * Mean(y));
}

// Correlation coefficient:
// (Pearson's product-moment coefficient)
template <class T> double Corr(vector<T>& x, vector<T>& y) {
  return Cov(x, y) / (Std(x) * Std(y));
}

// First requirment for weak stationary test, wild-sense stationary(WSS);
//  Expectation of {x_t} is a constant respect to t;
// VECTOR<VECTOR<DOUBLE> >& data: sampled data of {x_t},
//  data[i] stands for {x_i};
// VECTOR<DOUBLE> means: means[i] stands for the expectation of {x_i};
// Return: None;
void Means(vector<vector<double> > & data, vector<double>& means);

// Requirment for weak stationary test, wild-sense stationary(WSS);
//  Standard deviation of {x_t} is a constant respect to t;
// VECTOR<VECTOR<DOUBLE> >& data: sampled data of {x_t},
//  data[i] stands for {x_i};
// VECTOR<DOUBLE> means: means[i] stands for the expectation of {x_i};
// Return: None;
void Stds(vector<vector<double> > & data, vector<double>& stds);

// Second requirment for weak stationary test, wild-sense stationary(WSS);
//  covariance of {x_t} and {x_{t+\tau} is the function of \tau;
// VECTOR<VECTOR<DOUBLE> >& data: sampled data of {x_t},
//  data[i] stands for {x_i};
// VECTOR<DOUBLE> means: covs[i][j] stands for the between {x_j} and {x_{j + i + 1}};
// Return: None;
void Rule2(vector<vector<double> > & data, vector<vector<double> >& covs, size_t maxlag);

// Autocovariance of wide-sense stationary time series;
// VECTOR<double> data;
// VECTOR<double> ac:
// INT num_lag: number of timelag of autocovariance;
void AC(vector<double>& data, vector<double>& ac, size_t num_lag);

#endif
