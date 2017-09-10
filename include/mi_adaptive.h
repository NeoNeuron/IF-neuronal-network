//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-09-10
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

// Histogram of discrete bools;
// Return: the probability for TRUE;
double HistBool(vector<bool> & data);

// Joint probability function with adaptive partition of x and y; Iteration starts from the middle of the sorted sequence;
// VECTOR<DOUBLE> x: variable x;
// VECTOR<DOUBLE> y: variable y;
// VECTOR<DOUBLE> x_edges: edges for adaptive partition of x;
// VECTOR<DOUBLE> y_edges: edges for adaptive partition of y;
// VECTOR<VECTOR<DOUBLE> jointpdf: joint probability distribution function of x and y;
void JointPDF(vector<double>& x, vector<double>& y, vector<double>& x_edges, vector<double>& y_edges, vector<vector<double> >& jointpdf);

// Joint probability function with adaptive partition of x and y; Iteration starts from the middle of the sorted sequence;
// VECTOR<BOOL> binary_spikes: binary spike trains;
// VECTOR<DOUBLE> lfp: continuous local field potential;
// VECTOR<DOUBLE> lfp_edges: edges for adaptive partition of x;
// VECTOR<VECTOR<DOUBLE> jointpdf: joint probability distribution function of spike train and lfp;
void JointPDF(vector<bool>& binary_spikes, vector<double>& lfp, vector<double>& lfp_edges, vector<vector<double> >& jointpdf);

// Mutual information of two double sequences;
// VECTOR<DOUBLE> x, y: two original double sequences;
// Return: value of mutual information;
double MI(vector<double>& x, vector<double>& y);

// Mutual information between binary spiking train and LFP, with adaptive binning size (identical occupancy for each bin);
// VECTOR<BOOL> binary_spikes: original spike train with binary version;
// VECTOR<DOUBLE> LFP: local field potential;
// INT bin_number: number of bins of historgram of LFP;
// Return: valude of Mutual information;
double MI(vector<bool>& bool_series, vector<double>& double_series, int num_bin);


// Delayed mutual information of spiking train and LFP, with adaptive partitions;
void TDMI(vector<bool>& bool_series, vector<double>& double_series, int negative_time_delay, int positive_time_delay, vector<double> & tdmi, bool random_switch);

#endif // _IFNET_MI_H_
