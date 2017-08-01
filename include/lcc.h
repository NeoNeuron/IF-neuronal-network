#ifndef _IFNET_LCC_H_
#define _IFNET_LCC_H_

#include<string>
#include<vector>
using namespace std;

//	Read 1D data from *.txt files; Data type: (double);
//	Return: none;
// void ReadData(string path, vector<double> & data);

// Return mean value of elements of x;
template <class T> double Mean(vector<T>& x);

// Return standard deviation of elements of x;
template <class T> double Std(vector<T>& x);

// Pearson's product-moment coefficient: Linear correlated coefficient;
template <class T1, class T2> double LC(vector<T1>& x, vector<T2>& y);

void TDLC(vector<int>& raster, vector<double>& lfp, int negative_time_delay, int positive_time_delay, vector<double>& tdlc);

void TDLC(vector<double>& first, vector<double>& second, int negative_time_delay, int positive_time_delay, vector<double>& tdlc);

#endif // _IFNET_LCC_H_
