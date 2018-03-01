#include "../include/vecmanip.h"
#include <iostream>
using namespace std;

int main() {
  int n = 100;
  vector<int> shape(2);
  shape[0] = 10;
  shape[1] = 10;
  vector<int> data(n);
  for (size_t i = 0; i < data.size(); i ++) data[i] = i;
  for (vector<int>::iterator it = data.begin(); it != data.end(); it ++) {
    cout << *it << ',';
  }
  cout << '\n';
  vector<vector<int> > newdata;
  if (Reshape(data, newdata, shape)) {
    for (vector<vector<int> >::iterator it = newdata.begin(); it != newdata.end(); it ++) {
      for (vector<int>::iterator itt = it -> begin(); itt != it -> end(); itt ++) {
        cout << *itt << ',';
      }
      cout << '\n';
    }
  } else cout << "Failed\n";

  return 0;
}
