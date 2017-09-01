#ifndef _IFNET_CC_H_
#define _IFNET_CC_H_

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

// Pearson's product-moment coefficient: Calculate Pearson's correlation coefficient between two random variables, x and y;
double PCC(vector<double>& x, vector<double>& y);

// Calculate cross-correlation (or autocorrelation) of variables;
// VECTOR<DOUBLE> x: target (first) variable;
// VECTOR<VECTOR<DOUBLE> > y: reference (second) variable, number of delays by number of trials;
// VECTOR<DOUBLE> cc: container for cross-correlation or autocorrelation;
// Return: none;
void Correlation(vector<double>& x, vector<vector<double> >& y, vector<double>& cc);

#endif // _IFNET_CC_H_
