#include "../include/sta.h"
#include <iostream>
#include <numeric>

using namespace std;

double Product(vector<bool>& bool_series, vector<double>& double_series) {
  if (bool_series.size() != double_series.size()) {
    cout << "Series with inequal size." << endl;
    return 0.0;
  } else {
    int N = bool_series.size();
    vector<double> pdt(N, 0.0);
    for (int i = 0; i < N; i ++) {
      if (bool_series[i]) pdt[i] = double_series[i];
    }
    return accumulate(pdt.begin(), pdt.end(), 0.0);
  }
}
