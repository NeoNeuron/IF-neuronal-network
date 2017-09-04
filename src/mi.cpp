//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-02-27 21:48:11
//	Description: source file of mi.h
//***************
#include "../include/mi.h"
#include "../include/io.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <sstream>
using namespace std;

const double PI = 3.1415926;

bool comp(const double &a, const double &b) {
	return a < b;
}

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

// void int2str(const int &int_temp,string &string_temp) {
// 	stringstream stream;
// 	stream << int_temp;
// 	string_temp=stream.str();
// }

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

void JointPDF(vector<double>& x, vector<double>& y, double x_bin_width, double y_bin_width, vector<vector<double> >& jointpdf) {
	size_t num_pairs = x.size();
	double x_max = Max(x);
	double x_min = Min(x);
	double y_max = Max(y);
	double y_min = Min(y);
	size_t x_bin_number = ceil((x_max - x_min) / x_bin_width);
	size_t y_bin_number = ceil((y_max - y_min) / y_bin_width);
	jointpdf.resize(x_bin_number, vector<double>(y_bin_number, 0));
	int indx, indy;
	for (int i = 0; i < num_pairs; i++) {
		indx = floor((x[i] - x_min) / x_bin_width);
		indy = floor((y[i] - y_min) / y_bin_width);
		jointpdf[indx][indy] += 1.0 / num_pairs;
	}
}

void JointPDF(vector<bool>& binary_spikes, vector<double>& lfp, vector<double>& lfp_edges,
	vector<vector<double> >& jointpdf) {
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

double HistBool(vector<bool>& data) {
	int count = 0;
	for (vector<bool>::iterator it = data.begin(); it != data.end(); it++) {
		if (*it) count ++;
	}
	return count * 1.0 / data.size();
}

void HistInt(vector<int>& data, vector<double> & histogram, int min, int max) {
	int number = max - min + 1;
	vector<int> count(number, 0);
	int IND;
	for (vector<int>::iterator it = data.begin(); it != data.end(); it++) {
		IND = *it - min;
		count[IND] ++;
	}
	histogram.resize(number);
	for (int i = 0; i < number; i++) {
		histogram[i] = count[i] * 1.0 / data.size();
	}
}

void HistDouble(vector<double> & data, vector<double> & histogram, double bin_width) {
	double dmax = Max(data);
	double dmin = Min(data);
	size_t number = ceil((dmax -dmin) / bin_width);
	vector<int> count(number, 0);
	int IND;
	for (vector<double>::iterator it = data.begin(); it != data.end(); it++) {
		IND = floor((*it - dmin) / bin_width);
		if (IND == number) IND --;
		count[IND] ++;
	}
	histogram.resize(number);
	for (int i = 0; i < number; i++) {
		histogram[i] = count[i] * 1.0 / data.size();
	}
}

double MI(vector<bool>& x, vector<bool>& y) {
	if (x.size() != y.size()) {
		cout << "ERROR: x and y don't have the same length." << endl;
		return 0;
	} else {
		int np = x.size();
		// Histograms;
		double px, py;
		px = HistBool(x);
		py = HistBool(y);
		vector<double> Px(2), Py(2);
		Px[0] = 1 - px;
		Px[1] = px;
		Py[0] = 1 - py;
		Py[1] = py;

		vector<int> add(2, 0);
		vector<vector<int> > count_xy(2, add);
		for (int i = 0; i < np; i++) {
			if (x[i]) {
				if (y[i]) count_xy[1][1] ++;
				else count_xy[1][0] ++;
			} else {
				if (y[i]) count_xy[0][1] ++;
				else count_xy[0][0] ++;
			}
		}

		// Mutual information;
		double mi = 0;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				if (count_xy[i][j] != 0) {
					mi += count_xy[i][j] * 1.0 / np*log2(count_xy[i][j] * 1.0 / np / Px[i] / Py[j]);
				}
			}
		}
		return mi;
	}
}

// // Mutual information with uniformly binning histogram;
// double MI(vector<double>& x, vector<double>& y, double x_bin_width, double y_bin_width) {
// 	//	Compare the length of x & y;
// 	int length;
// 	vector<double> x_copy = x, y_copy = y;
// 	if (x_copy.size() > y_copy.size()) {
// 		length = y_copy.size();
// 		x_copy.erase(x_copy.begin() + length, x_copy.end());
// 	} else {
// 		length = x_copy.size();
// 		y_copy.erase(y_copy.begin() + length, y_copy.end());
// 	}
// 	// Calculate actual occupancy;
// 	// Find edges of histogram;
// 	vector<double> x_hist, y_hist;
// 	HistDouble(x_copy, x_hist, x_bin_width);
// 	HistDouble(y_copy, y_hist, y_bin_width);
//
// 	// cout << x_hist.size() << ',' << y_hist.size() << '\n';
//
// 	// 	Calculate conditional probability;
// 	vector<vector<double> > jointpdf;
// 	JointPDF(x_copy, y_copy, bin_width, bin_width, jointpdf);
//
// 	// Calculate mutual information;
// 	double mi = 0;
// 	for (size_t i = 0; i < x_hist.size(); i++) {
// 		for (size_t j = 0; j < y_hist.size(); j++) {
// 			if (jointpdf[i][j] != 0) {
// 				mi += jointpdf[i][j] * log2(jointpdf[i][j] / (x_hist[i] * y_hist[j]));
// 			}
// 		}
// 	}
// 	return mi;
// }

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
	} else {
		length = x.size();
	}
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
	double mi = 0;
	for (int i = 0; i < bin_number; i++) {
		for (int j = 0; j < bin_number; j++) {
			if (jointpdf[i][j] != 0) {
				mi += jointpdf[i][j] * log(pow(bin_number, 2)*jointpdf[i][j]);
			}
		}
	}
	return mi;
}

double MI(vector<bool>& bool_series, vector<double>& double_series, int num_bin) {
	// calculate the occupancy;
	int occupancy = double_series.size() / num_bin;
	int num_pair = occupancy * num_bin;
	vector<double> double_copy = double_series;
	double_copy.erase(double_copy.begin() + num_pair, double_copy.end());

	// calculate histogram;
	vector<double> edges;
	FindEdges(double_copy, edges, occupancy);
	double p_spike;
	p_spike = HistBool(bool_series);

	// Calculate conditional probability;
	vector<vector<double> > jointpdf;
	JointPDF(bool_series, double_copy, edges, jointpdf);

	// Calculate mutual information;
	double mi = 0; // Pxy = P(x,y);
	for (int i = 0; i< num_bin; i++) {
		if (jointpdf[0][i] != 0) {
			mi += jointpdf[0][i] *log2(num_bin*jointpdf[0][i] / p_spike);
		}
	}
	for (int i = 0; i< num_bin; i++) {
		if (jointpdf[1][i] != 0) {
			mi += jointpdf[1][i] *log2(num_bin*jointpdf[1][i] / (1 - p_spike));
		}
	}
	return mi;
}

void TDMI(vector<bool>& x, vector<bool>& y, double dt, double tmax, int negative_time_delay, int positive_time_delay, vector<double> & tdmi) {
	tdmi.resize(positive_time_delay + negative_time_delay + 1, 0);

	// No shift;
	tdmi[negative_time_delay] = MI(x, y);

	// Negative shift;
	vector<bool> x_copy = x, y_copy = y;
	for (int i = 0; i < negative_time_delay; i++) {
		x_copy.erase(x_copy.begin());
		y_copy.erase(y_copy.end() - 1);
		tdmi[negative_time_delay - i - 1] = MI(x_copy, y_copy);
	}

	// Positive shift;
	x_copy = x;
	y_copy = y;
	for (int i = 0; i < positive_time_delay; i++) {
		x_copy.erase(x_copy.end() - 1);
		y_copy.erase(y_copy.begin());
		tdmi[negative_time_delay + i + 1] = MI(x_copy, y_copy);
	}
}

void TDMI(vector<double>& x, vector<vector<double> >& y, vector<double> & tdmi) {
	tdmi.resize(y.size(), 0);
	for (size_t i = 0; i < y.size(); i++) {
		tdmi[i] = MI(x, y[i]);
	}
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
