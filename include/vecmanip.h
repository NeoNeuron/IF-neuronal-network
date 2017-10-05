#ifndef _IFNET_VECMANIP_H_
#define _IFNET_VECMANIP_H_
#include <vector>
#include <string>
using namespace std;

// Transpose 2D vector;
// VECTOR<T> data: original data;
// VECTOR<T> newdata: transposed data;
template <class T> bool Transpose(vector<T>& data, vector<T>& newdata);

// Muliple all elements in vectors;
// VECTOR<T> data: original data;
// Return: value of result of Pi;
template <class T> T Pi(vector<T> data);

// Reshape 1D vector to multi-dimensional vectors;
// VECTOR<T> data: original data;
// VECTOR<T> newdata: new data;
// SIZE_T* shape: sizes of 2-d vectors;
template <class T> bool Reshape(vector<T>& data, vector<T>& newdata, size_t* shape);

#include "vecmanip.hpp"

#endif
