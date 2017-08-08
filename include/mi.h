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

//	Gaussian kernel:
//	Return: random number obey the standard Gaussion distribution;
double GaussKernel();

//	Poisson generator:
//	DOUBLE rate: mean Poisson rate;
//	DOUBLE t_last: last time point for Poisson events;
//	Return: time point for next Poisson events;
double PoissonGenerator(double rate, double t_last);

//	Find Edges: Find edges of bins with adaptive partition, ie. each bin is equally occupied;
//	VECTOR<DOUBLE> data: original data;
//	VECTOR<DOUBLE> edges: edges of bins of histogram;
//	INT occupancy: number of data points occupied in each bin;
//	INT residue: number of data points at the end of the vector which is not considered into histogram;
//	Return: none;
void FindEdges(vector<double>& data, vector<double>& edges, int occupancy);

//  Joint probability function with adaptive partition of x and y; Iteration starts from the middle of the sorted sequence;
//  VECTOR<DOUBLE> x: variable x;
//  VECTOR<DOUBLE> y: variable y;
//  VECTOR<DOUBLE> x_edges: edges for adaptive partition of x;
//  VECTOR<DOUBLE> y_edges: edges for adaptive partition of y;
//  VECTOR<VECTOR<DOUBLE> jointpdf: joint probability distribution function of x and y;
void JointPDF(vector<double>& x, vector<double>& y, vector<double>& x_edges, vector<double>& y_edges, vector<vector<double> >& jointpdf);

//  Joint probability function with adaptive partition of x and y; Iteration starts from the middle of the sorted sequence;
//  VECTOR<BOOL> binary_spikes: binary spike trains;
//  VECTOR<DOUBLE> lfp: continuous local field potential;
//  VECTOR<DOUBLE> lfp_edges: edges for adaptive partition of x;
//  VECTOR<VECTOR<DOUBLE> jointpdf: joint probability distribution function of spike train and lfp;
void JointPDF(vector<bool>& binary_spikes, vector<double>& lfp, vector<double>& lfp_edges, vector<vector<double> >& jointpdf);

// Find maximum value in data;
// Return max;
template <class T> T Max(vector<T>& data);

// Find minimum value in data;
// Return min;
template <class T> T Min(vector<T>& data);

//	Convert double spike train to binary sequence;
//	VECTOR<DOUBLE> spikes: original spike train;
//	VECTOR<BOOL> binary_spikes: binary sequence; true for having spike; false for no spike;
//	DOUBLE tmax: maximum time of time range;
//	DOUBLE dt: the size of time step that a binary value considered;
//	Return: none;
void Spike2Bool(vector<double> & spikes, vector<bool> & binary_spikes, double tmax, double dt);

// 	Histogram of discrete bools;
//	Return: the probability for TRUE;
double HistBool(vector<bool> & data);

// 	Histogram of discrete ints;
//	VECTOR<INT> data: original int data;
//	VECTOR<DOUBLE> histogram: histogram of int;
//	INT min: lower limit of int data; [minimum value]
//	INT max: upper limit of int data; [maximum value]
//	Return: none;
void HistInt(vector<int> & data, vector<double> & histogram, int min, int max);

// 	Histogram of continues variables;
//	VECTOR<DOUBLE> data: original double data;
//	VECTOR<DOUBLE> histogram: histogram of double;
//	DOUBLE min: lower limit of int data; [minimum value]
//	DOUBLE max: upper limit of int data; [maximum value]
//	DOUBLE bin_width: bin width of histogram; it is uniformly binned;
//	Return: none;
void HistDouble(vector<double> & data, vector<double> & histogram, double min, double max, double bin_width);

// 	Mutual information of two binary spike trains;
//	VECTOR<BOOL> x, y: tow original bool sequence;
//	Return: value of mutual information;
double MI(vector<bool>& x, vector<bool>& y);

// 	Mutual information of two double sequences;
//	VECTOR<DOUBLE> x, y: two original double sequences;
//	DOUBLE* x_max_and_min, y_max_and_min: set of maximum and minimum value; *_max_and_min[0] = *_max; *_max_and_min[1] = *_min;
//	DOUBLE x_bin_width, y_bin_width: bin width of histograms of x and y;
//	Return: value of mutual information;
double MI(vector<double>& x, vector<double>& y, double* x_max_and_min, double* y_max_and_min, double x_bin_width, double y_bin_width);

// 	Mutual information of two double sequences;
//	VECTOR<DOUBLE> x, y: two original double sequences;
//	Return: value of mutual information;
double MI(vector<double>& x, vector<double>& y);
double MI(vector<double>& x, vector<double>& y, string pdf_path);

// 	Mutual information between binary spiking train and LFP; [Based on histogram scheme which has identical occupancy for each bin;
//	VECTOR<BOOL> binary_spikes: original spike train with binary version;
//	VECTOR<DOUBLE> LFP: local field potential;
//	INT bin_number: number of bins of historgram of LFP;
//	Return: valude of Mutual information;
double MI(vector<bool> &binary_spikes, vector<double> &LFP, int bin_number);

// Time-delayed mutual information between two spike trains;
void TDMI(vector<double>& x, vector<double>& y, double dt, double tmax, int negative_time_delay, int positive_time_delay, vector<double> & tdmi);

// Delayed mutual information of two double vector, by applying histogram scheme with equal bin size;
void TDMI_uniform(vector<double>& x, vector<double>& y, int negative_time_delay, int positive_time_delay, vector<double> & tdmi);

// Delayed mutual information of two double vector, wiht adaptive partition;
void TDMI_adaptive(vector<double>& x, vector<double>& y, int negative_time_delay, int positive_time_delay, vector<double> & tdmi);

// Delayed mutual information of spiking train and LFP, with adaptive partitions;
void TDMI(vector<double>& spikes, vector<double>& LFP, double dt, double sampling_dt, int negative_time_delay, int positive_time_delay, vector<double> & tdmi, bool random_switch);

#endif // _IFNET_MI_H_
