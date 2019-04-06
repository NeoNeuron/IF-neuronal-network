//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-09-10
//	Description: source file of mi_uniform.h
//***************
#include "../include/mi_adaptive.h"
#include <algorithm>
#include <cmath>
using namespace std;

bool comp(const double &a, const double &b) { return a < b; }

double HistBool(vector<bool>& data) {
	int num_true = count(data.begin(), data.end(), true);
	return num_true * 1.0 / data.size();
}

void FindEdges(vector<double>& data, vector<double>& edges, int occupancy) {
	vector<double> data_copy = data;
	sort(data_copy.begin(), data_copy.end(), comp);
	size_t N = floor(data_copy.size() / occupancy);
	edges.resize(N);
	for (size_t i = 0; i < N; i++) {
		edges[i] = data_copy[i*occupancy];
	}
}

void JointPDF(vector<double>& x, vector<double>& y, vector<double>& x_edges, vector<double>& y_edges, vector<vector<double> >& jointpdf) {
	int num_pairs = x.size();
	int x_bin_number = x_edges.size();
	int y_bin_number = y_edges.size();
	jointpdf.resize(x_bin_number, vector<double>(y_bin_number, 0));
	int indx, indy;
	for (int i = 0; i < num_pairs; i++) {
		indx = floor(x_bin_number / 2);
		if (x[i] < x_edges[indx]) {
			for (int j = indx - 1; j > -1; j--) {
				if (x[i] >= x_edges[j]) {indx = j; break;}
			}
		} else {
			for (int j = indx + 1; j < x_bin_number; j++) {
				if (x[i] < x_edges[j]) {indx = j - 1; break;}
			}
			if (indx == floor(x_bin_number / 2)) indx = x_bin_number - 1;
		}
		indy = floor(y_bin_number / 2);
		if (y[i] < y_edges[indy]) {
			for (int j = indy - 1; j > -1; j--) {
				if (y[i] >= y_edges[j]) {indy = j; break;}
			}
		} else {
			for (int j = indy + 1; j < y_bin_number; j++) {
				if (y[i] < y_edges[j]) {indy = j - 1; break;}
			}
			if (indy == floor(y_bin_number / 2)) indy = y_bin_number - 1;
		}
		jointpdf[indx][indy] ++;
	}
	for (vector<vector<double> >::iterator it = jointpdf.begin(); it != jointpdf.end(); it++) {
		for (vector<double>::iterator itt = it->begin(); itt != it->end(); itt++) {
			*itt /= num_pairs;
		}
	}
}

void JointPDF(vector<bool>& binary_spikes, vector<double>& lfp, vector<double>& lfp_edges, vector<vector<double> >& jointpdf) {
	int bin_number = lfp_edges.size();
	int num_pairs = lfp.size();
	jointpdf.resize(2, vector<double>(bin_number, 0));
	int ind;
	for (int i = 0; i < num_pairs; i++) {
		// determine lfp's coordination;
		ind = floor(bin_number / 2);
		if (lfp[i] < lfp_edges[ind]) {
			for (int j = ind - 1; j > -1; j--) {
				if (lfp[i] >= lfp_edges[j]) {ind = j; break;}
			}
		} else {
			for (int j = ind + 1; j < bin_number; j++) {
				if (lfp[i] < lfp_edges[j]) {ind = j - 1; break;}
			}
			if (ind == floor(bin_number / 2)) ind = bin_number - 1;
		}
		// determine spike's coordination;
		if (binary_spikes[i]) jointpdf[0][ind]++;
		else jointpdf[1][ind]++;
	}
	for (vector<vector<double> >::iterator it = jointpdf.begin(); it != jointpdf.end(); it++) {
		for (vector<double>::iterator itt = it->begin(); itt != it->end(); itt++) {
			*itt /= num_pairs;
		}
	}
}

double MI(vector<double>& x, vector<double>& y) {
	//	Compare the length of x & y;
	int length;
	vector<double> x_copy = x, y_copy = y;
	if (x.size() > y.size()) {
		length = y.size();
		x_copy.erase(x_copy.begin() + length, x_copy.end());
	} else if (x.size() < y.size()) {
		length = x.size();
		y_copy.erase(y_copy.begin() + length, y_copy.end());
	} else length = x.size();
	// Calculate actual occupancy;
	int bin_number = floor(sqrt(length / 5)); // bin_number: numberber of non-uniform bins;
	int occupancy = floor(length / bin_number);
	// Find edges of histogram;
	vector<double> x_edges, y_edges;
	FindEdges(x_copy, x_edges, occupancy);
	FindEdges(y_copy, y_edges, occupancy);
	// 	Calculate conditional probability;
	vector<vector<double> > jointpdf;
	JointPDF(x_copy, y_copy, x_edges, y_edges, jointpdf);
	// Calculate mutual information;
	double mi = 0.0;
	for (int i = 0; i < bin_number; i++) {
		for (int j = 0; j < bin_number; j++) {
			if (abs(jointpdf[i][j]) < 1e-20) {
				mi += jointpdf[i][j] * log(pow(bin_number, 2)*jointpdf[i][j]);
			}
		}
	}
	return mi;
}

double MI(vector<bool>& bool_series, vector<double>& double_series, int bin_num) {
	// calculate the occupancy;
	int occupancy = double_series.size() / bin_num;
	int num_pair = occupancy * bin_num;
	vector<double> double_copy = double_series;
	double_copy.erase(double_copy.begin() + num_pair, double_copy.end());
	// calculate histogram;
	vector<double> edges;
	FindEdges(double_copy, edges, occupancy);
	double p_spike[2];
	p_spike[0] = HistBool(bool_series);
	p_spike[1] = 1 - p_spike[0];
	// Calculate conditional probability;
	vector<vector<double> > jointpdf;
	JointPDF(bool_series, double_copy, edges, jointpdf);
	// Calculate mutual information;
	double mi = 0.0; // Pxy = P(x,y);
	for (int i = 0; i< bin_num; i++) {
		for (int j = 0; j < bin_num; j++) {
			if (abs(jointpdf[i][j]) < 2e-20) {
				mi += jointpdf[i][j] *log(bin_num*jointpdf[i][j] / p_spike[i]);
			}
		}
	}
	return mi;
}

void TDMI(vector<bool>& bool_series, vector<double>& double_series, int negative_time_delay, int positive_time_delay, vector<double>& tdmi, bool random_switch) {
	int num_reps = negative_time_delay + positive_time_delay + 1;
	tdmi.resize(num_reps, 0);
	int N = bool_series.size(); // number of reduced LFP data point;
	int num_spike = count(bool_series.begin(), bool_series.end(), true);
	int num_bin;
	if (num_spike != 0) {
		if (num_spike < 5) num_bin = 1;
		else num_bin = floor(num_spike / 5); // expected_occupancy = 5;
		if (random_switch == true) {
			random_shuffle(bool_series.begin(), bool_series.end());
		}
		// No shift;
		tdmi[negative_time_delay] = MI(bool_series, double_series, num_bin);

		// Negative shift;
		vector<bool> bool_copy = bool_series;
		vector<double> double_copy = double_series;
		for (int i = 0; i < negative_time_delay; i++) {
			bool_copy.erase(bool_copy.begin());
			double_copy.erase(double_copy.end() - 1);
			tdmi[negative_time_delay - i - 1] = MI(bool_copy, double_copy, num_bin);
		}

		// Positive shift;
		bool_copy = bool_series;
		double_copy = double_series;
		for (int i = 0; i < positive_time_delay; i++) {
			bool_copy.erase(bool_copy.end() - 1);
			double_copy.erase(double_copy.begin());
			tdmi[negative_time_delay + i + 1] = MI(bool_copy, double_copy, num_bin);
		}
	}
}
