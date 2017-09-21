#include "vecmanip.h"
#include <iostream>
using namespace std;

template <class T> bool Transpose(vector<vector<T> >& data, vector<vector<T> >& newdata) {
  // Judge the transpose requirement;
  size_t firstlen = data.begin()->size();
  bool mat_set = true;
  for (typename vector<vector<T> >::iterator it = data.begin() + 1; it != data.end(); it ++) {
    if (it->size() != firstlen) {
      mat_set = false;
      break;
    }
  }
  if (mat_set) {
    newdata.resize(firstlen, vector<T>(data.size()));
    for (size_t i = 0; i < data.size(); i ++){
      for (size_t j = 0; j < firstlen; j ++) {
        newdata[j][i] = data[i][j];
      }
    }
  }
  return mat_set;
}

template <class T> T Pi(vector<T> data) {
  T pi = 1;
  for (typename vector<T>::iterator it = data.begin(); it != data.end(); it++) {
    pi = pi * (*it);
  }
  return pi;
}

template <class T> bool Reshape(vector<T>& data, vector<vector<T> >& newdata, vector<int>& shape) {
  size_t length = Pi(shape);
  if (data.size() != length) return false;
  else if (shape.size() == 2) {
    newdata.resize(shape[0], vector<T>(shape[1]));
    for (size_t i = 0; i < shape[0]; i ++) {
      for (size_t j = 0; j < shape[1]; j ++) {
        newdata[i][j] = data[i*shape[0] + j];
      }
    }
    return true;
  } else return false;
}
