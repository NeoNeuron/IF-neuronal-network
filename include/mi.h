//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-02-27 21:48:02
//	Description: Define algorithm to estimate Mutual Information(MI) and Time-Delayed Mutual Information(TDMI);
//***************
#ifndef _IFNET_MI_H_
#define _IFNET_MI_H_

#include <string>
#include <vector>
using namespace std;

// Find Edges: Find edges of bins with adaptive partition, ie. each bin is equally occupied;
// VECTOR<DOUBLE> data: original data;
// VECTOR<DOUBLE> edges: edges of bins of histogram;
// INT occupancy: number of data points occupied in each bin;
// INT residue: number of data points at the end of the vector which is not considered into histogram;
// Return: none;
void FindEdges(vector<double>& data, vector<double>& edges, int occupancy);

// Joint probability function with adaptive partition of x and y; Iteration starts from the middle of the sorted sequence;
// VECTOR<DOUBLE> x: variable x;
// VECTOR<DOUBLE> y: variable y;
// VECTOR<DOUBLE> x_edges: edges for adaptive partition of x;
// VECTOR<DOUBLE> y_edges: edges for adaptive partition of y;
// VECTOR<VECTOR<DOUBLE> jointpdf: joint probability distribution function of x and y;
void JointPDF(vector<double>& x, vector<double>& y, vector<double>& x_edges, vector<double>& y_edges, vector<vector<double> >& jointpdf);

void JointPDF(vector<bool>& x, vector<bool>& y, vector<vector<double> >& jointpdf);

void JointPDF(vector<double>& x, vector<double>& y, size_t x_bin_num, size_t y_bin_num, vector<vector<double> >& jointpdf);

// Joint probability function with uniform partition of x and y;
// VECTOR<BOOL> binary_spikes: binary spike trains;
// VECTOR<DOUBLE> lfp: continuous local field potential;
// DOUBLE bin_width: width of uniform partition of x;
// VECTOR<VECTOR<DOUBLE> jointpdf: joint probability distribution function of spike train and lfp;
void JointPDF(vector<bool>& binary_spikes, vector<double>& lfp, double bin_num, vector<vector<double> >& jointpdf);

// Joint probability function with adaptive partition of x and y; Iteration starts from the middle of the sorted sequence;
// VECTOR<BOOL> binary_spikes: binary spike trains;
// VECTOR<DOUBLE> lfp: continuous local field potential;
// VECTOR<DOUBLE> lfp_edges: edges for adaptive partition of x;
// VECTOR<VECTOR<DOUBLE> jointpdf: joint probability distribution function of spike train and lfp;
void JointPDF(vector<bool>& binary_spikes, vector<double>& lfp, vector<double>& lfp_edges, vector<vector<double> >& jointpdf);

// Marginal probability density function of bivariable joint distribution density function;
// Return: none;
void MarginalPDF(vector<vector<double> >& jointpdf, vector<double>& p1, vector<double>& p2);

// Histogram of discrete bools;
// Return: the probability for TRUE;
double HistBool(vector<bool> & data);

// Histogram of discrete ints;
// VECTOR<INT> data: original int data;
// VECTOR<DOUBLE> histogram: histogram of int;
// INT min: lower limit of int data; [minimum value]
// INT max: upper limit of int data; [maximum value]
// Return: none;
void HistInt(vector<int> & data, vector<double> & histogram, int min, int max);

// Histogram of continues variables;
// VECTOR<DOUBLE> data: original double data;
// VECTOR<DOUBLE> histogram: histogram of double;
// DOUBLE min: lower limit of int data; [minimum value]
// DOUBLE max: upper limit of int data; [maximum value]
// DOUBLE bin_width: bin width of histogram; it is uniformly binned;
// Return: none;
void HistDouble(vector<double> & data, vector<double> & histogram, double bin_width);

// Mutual information of two binary spike trains;
// VECTOR<BOOL> x, y: tow original bool sequence;
// Return: value of mutual information;
double MI(vector<bool>& x, vector<bool>& y);

// Mutual information of two double sequences;
// VECTOR<DOUBLE> x, y: two original double sequences;
// DOUBLE* x_max_and_min, y_max_and_min: set of maximum and minimum value; *_max_and_min[0] = *_max; *_max_and_min[1] = *_min;
// DOUBLE x_bin_width, y_bin_width: bin width of histograms of x and y;
// Return: value of mutual information;
double MI(vector<double>& x, vector<double>& y, double x_bin_width, double y_bin_width);

// Mutual information of two double sequences;
// VECTOR<DOUBLE> x, y: two original double sequences;
// Return: value of mutual information;
double MI(vector<double>& x, vector<double>& y);

// Mutual information between binary spiking train and LFP with uniform binning size;
// VECTOR<BOOL> binary_spikes: original spike train with binary version;
// VECTOR<DOUBLE> LFP: local field potential;
// INT bin_number: number of bins of historgram of LFP;
// Return: valude of Mutual information;
double MI(vector<bool>& bool_series, vector<double>& double_series, int num_bin);

// // Mutual information between binary spiking train and LFP, with adaptive binning size (identical occupancy for each bin);
// // VECTOR<BOOL> binary_spikes: original spike train with binary version;
// // VECTOR<DOUBLE> LFP: local field potential;
// // INT bin_number: number of bins of historgram of LFP;
// // Return: valude of Mutual information;
// double MI(vector<bool>& bool_series, vector<double>& double_series, int num_bin);

// Time-delayed mutual information between two spike trains;
void TDMI(vector<double>& x, vector<double>& y, double dt, double tmax, int negative_time_delay, int positive_time_delay, vector<double> & tdmi);

// Delayed mutual information of two double series, wiht adaptive partition;
void TDMI(vector<double>& x, vector<vector<double> >& y, vector<double>& tdmi);

// Delayed mutual information of spiking train and LFP, with adaptive partitions;
void TDMI(vector<bool>& bool_series, vector<vector<double> >& double_series, vector<double> & tdmi, size_t bin_num);

// // Delayed mutual information of spiking train and LFP, with adaptive partitions;
// void TDMI(vector<bool>& bool_series, vector<double>& double_series, int negative_time_delay, int positive_time_delay, vector<double> & tdmi, bool random_switch);

#endif // _IFNET_MI_H_
