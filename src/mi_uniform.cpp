//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-09-10
//	Description: source file of mi_uniform.h
//***************
#include "../include/mi_uniform.h"
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
using namespace std;

double Max(vector<double>& data) {
	vector<double>::iterator it;
	it = max_element(data.begin(), data.end());
	return *it;
}

double Min(vector<double>& data) {
	vector<double>::iterator it;
	it = min_element(data.begin(), data.end());
	return *it;
}

void JointPDF(vector<bool>& x, vector<bool>& y, vector<vector<double> >& jointpdf) {
	size_t num_pairs = x.size();
	vector<vector<int> > count_xy(2, vector<int>(2, 0));
	for (int i = 0; i < num_pairs; i++) {
		if (x[i]) {
			if (y[i]) count_xy[1][1] ++;
			else count_xy[1][0] ++;
		} else {
			if (y[i]) count_xy[0][1] ++;
			else count_xy[0][0] ++;
		}
	}
	jointpdf.resize(2, vector<double>(2, 0.0));
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
		jointpdf[i][j] +=  count_xy[i][j] * 1.0 / num_pairs;
		}
	}
}

void JointPDF(vector<double>& x, vector<double>& y, size_t x_bin_num, size_t y_bin_num, vector<vector<double> >& jointpdf) {
	size_t num_pairs = x.size();
	double x_max = Max(x);
	double x_min = Min(x);
	double y_max = Max(y);
	double y_min = Min(y);
	double x_bin_width = (x_max - x_min) / x_bin_num;
	double y_bin_width = (y_max - y_min) / y_bin_num;
	jointpdf.resize(x_bin_num, vector<double>(y_bin_num, 0.0));
	int indx, indy;
	for (int i = 0; i < num_pairs; i++) {
		indx = floor((x[i] - x_min) / x_bin_width);
		indy = floor((y[i] - y_min) / y_bin_width);
		if (indx == x_bin_num) indx = x_bin_num - 1;
		if (indy == y_bin_num) indy = y_bin_num - 1;
		jointpdf[indx][indy] += 1.0 / num_pairs;
	}
}

void JointPDF(vector<bool>& binary_spikes, vector<double>& lfp, size_t bin_num, vector<vector<double> >& jointpdf) {
	size_t num_pairs = lfp.size();
	jointpdf.resize(2, vector<double>(bin_num, 0.0));
	double lfp_min = Min(lfp);
	double bin_width = (Max(lfp) - lfp_min) / bin_num;
	int ind;
	for (size_t i = 0; i < num_pairs; i++) {
		// determine lfp's coordination;
		ind = floor((lfp[i] - lfp_min) / bin_width);
		if (ind == bin_num) ind = bin_num - 1;
		// determine spike's coordination;
		if (binary_spikes[i]) jointpdf[1][ind] += 1.0 / num_pairs;
		else jointpdf[0][ind] += 1.0 / num_pairs;
	}
}

void MarginalPDF(vector<vector<double> >& jointpdf, vector<double>& p1, vector<double>& p2) {
	p1.resize(jointpdf.size());
	p2.clear();
	p2.resize(jointpdf.begin() -> size(), 0.0);
	for (size_t i = 0; i < p1.size(); i++) {
		p1[i] = accumulate(jointpdf[i].begin(), jointpdf[i].end(), 0.0);
		for (size_t j = 0; j < p2.size(); j++) {
			p2[j] += jointpdf[i][j];
		}
	}
}

double MI(vector<bool>& x, vector<bool>& y) {
	if (x.size() != y.size()) {
		cout << "ERROR: x and y don't have the same length." << endl;
		return 0;
	} else {
		// Joint Probability Histograms;
		vector<vector<double> > joint_xy;
		JointPDF(x, y, joint_xy);
		// Marginal distribution;
		vector<double> px, py;
		MarginalPDF(joint_xy, px, py);
		// Mutual information;
		double mi = 0.0;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				if (abs(joint_xy[i][j]) < 1e-20) {
					mi += joint_xy[i][j]*log(joint_xy[i][j] / (px[i] * py[j]));
				}
			}
		}
		return mi;
	}
}

// Mutual information with uniformly binning histogram;
double MI(vector<double>& x, vector<double>& y, size_t x_bin_num, size_t y_bin_num) {
	//	Compare the length of x & y;
	if (x.size() != y.size()) {
		cout << "ERROR: x and y don't have the same length." << endl;
		return 0;
	} else {
		// Calculate joint probability distribution function;
		vector<vector<double> > jointpdf;
		JointPDF(x, y, x_bin_num, y_bin_num, jointpdf);
		// Marginal probability distribution function;
		vector<double> px, py;
		MarginalPDF(jointpdf, px, py);
		// Calculate mutual information;
		double mi = 0.0;
		for (size_t i = 0; i < px.size(); i++) {
			for (size_t j = 0; j < py.size(); j++) {
				if (abs(jointpdf[i][j]) < 1e-20) {
					mi += jointpdf[i][j] * log(jointpdf[i][j] / (px[i] * py[j]));
				}
			}
		}
		if (abs(mi) < 1e-15) mi = 0;
		return mi;
	}
}

double MI(vector<bool>& bool_series, vector<double>& double_series, int bin_num) {
	if (bool_series.size() != double_series.size()) {
		cout << "ERROR: x and y don't have the same length." << endl;
		return 0;
	} else {
		// calculate the occupancy;
		size_t num_pair = double_series.size();
		double bin_width = (Max(double_series) - Min(double_series)) / bin_num;
		// Calculate joint probability distribution function;
		vector<vector<double> > jointpdf;
		JointPDF(bool_series, double_series, bin_num, jointpdf);
		// calculate histogram;
		vector<double> p_spike, p_double;
		MarginalPDF(jointpdf, p_spike, p_double);
		// Calculate mutual information;
		double mi = 0.0; // Pxy = P(x,y);
		for (int i = 0; i< 2; i++) {
			for (int j = 0; j< bin_num; j++) {
				if (jointpdf[i][j] != 0) {
					mi += jointpdf[i][j] *log(jointpdf[i][j] / (p_spike[i] * p_double[j]));
				}
			}
		}
		if (abs(mi) < 1e-15) mi = 0;
		return mi;
	}
}

void TDMI(vector<bool>& x, vector<vector<bool> >& y, vector<double> & tdmi) {
	tdmi.resize(y.size(), 0);
	for (size_t i = 0; i < y.size(); i++) {
		tdmi[i] = MI(x, y[i]);
	}
}

void TDMI(vector<double>& x, vector<vector<double> >& y, vector<double> & tdmi, size_t x_bin_num, size_t y_bin_num) {
	tdmi.resize(y.size(), 0);
	for (size_t i = 0; i < y.size(); i++) {
		tdmi[i] = MI(x, y[i], x_bin_num, y_bin_num);
	}
}

void TDMI(vector<bool>& bool_series, vector<vector<double> >& double_series, vector<double>& tdmi, size_t bin_num) {
	size_t reps_num = double_series.size() + 1;
	tdmi.resize(reps_num, 0);
	for (size_t i = 0; i < double_series.size(); i++) {
		tdmi[i] = MI(bool_series, double_series[i], bin_num);
	}
}
