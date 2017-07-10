#ifndef _IFNET_LCC_H_
#define _IFNET_LCC_H_

#include<string>
#include<vector>
using namespace std;

//	Read 1D data from *.txt files; Data type: (double);
//	Return: none;
// void ReadData(string path, vector<double> & data);

// Return mean value of elements of x;
double Mean(vector<int>& x);
double Mean(vector<double>& x);

// Return standard deviation of elements of x;
double Std(vector<int>& x);
double Std(vector<double>& x);

// Pearson's product-moment coefficient: Linear correlated coefficient;
double LC(vector<int>& raster, vector<double>& lfp);

double LC(vector<double>& first, vector<double>& second);

void TDLC(vector<int>& raster, vector<double>& lfp, int negative_time_delay, int positive_time_delay, vector<double>& tdlc);

void TDLC(vector<double>& first, vector<double>& second, int negative_time_delay, int positive_time_delay, vector<double>& tdlc);

#endif // _IFNET_LCC_H_
