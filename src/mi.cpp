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

void int2str(const int &int_temp,string &string_temp) {
	stringstream stream;
	stream << int_temp;
	string_temp=stream.str();
}

double GaussKernel() {
	double U = -log(1 - rand() / (RAND_MAX*1.0));
	double V = 2 * PI * rand() / (RAND_MAX*1.0);
	double number = sqrt(2 * U)*cos(V);
	while (abs(number) > 1e10) {
		U = -log(1 - rand() / (RAND_MAX*1.0));
		V = 2 * PI * rand() / (RAND_MAX*1.0);
		number = sqrt(2 * U)*cos(V);
	}
	return number;
}

double PoissonGenerator(double rate, double t_last) {
	double x;
	x = rand() / (RAND_MAX + 1.0);
	while (abs(x) > 1e10) x = rand() / (RAND_MAX + 1.0);
	return t_last - log(x) / rate;
}

void FindEdges(vector<double>& data, vector<double>& edges, int occupancy, int residue) {
	vector<double> data_copy = data;
	data_copy.erase(data_copy.end() - residue, data_copy.end());
	sort(data_copy.begin(), data_copy.end(), comp);
	int N = floor(data_copy.size() / occupancy);
	edges.resize(N);
	for (int i = 0; i < N; i++) {
		edges[i] = data_copy[i*occupancy];
	}
	data_copy.clear();
}

void JointPDF(vector<double>& x, vector<double>& y, vector<double>& x_edges, vector<double>& y_edges, vector<vector<double> >& jointpdf) {
	int num_pairs = x.size();
	int x_bin_number = x_edges.size();
	int y_bin_number = y_edges.size();
	vector<double> add(y_bin_number, 0);
	jointpdf.resize(x_bin_number, add);
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

void JointPDF(vector<bool>& binary_spikes, vector<double>& lfp, vector<double>& lfp_edges,
	vector<vector<double> >& jointpdf) {
	int bin_number = lfp_edges.size();
	int num_pairs = lfp.size();
	jointpdf.resize(2);
	for (int i = 0; i < 2; i++) jointpdf[i].resize(lfp_edges.size(), 0);
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

void FindMaxMin(vector<double>& data, double *max_and_min, double bin_width) {
	double min, max;
	min = Min(data);
	max = Max(data);
	min = floor(min / bin_width)*bin_width;
	max = ceil(max / bin_width)*bin_width;
	*max_and_min = max;
	*(max_and_min + 1) = min;
}

void ConvertSpikeToBinary(vector<double>& spikes, vector<bool> & binary_spikes, double tmax, double dt) {
	int T = ceil(tmax / dt);
	binary_spikes.resize(T, false);
	for (vector<double>::iterator it = spikes.begin(); it != spikes.end(); it++) {
		int index = floor(*it / dt);
		if (index < T - 1) {
			binary_spikes[index] = true;
		}
	}
}

void ConvertBinaryToInt(vector<bool> & binary_spikes, vector<int> & spikeInt, int bins) {
	spikeInt.clear();
	int T = floor(binary_spikes.size() / bins);
	int add_int;
	for (int i = 0; i < T; i++) {
		add_int = 0;
		for (int j = 0; j < bins; j++) {
			if (binary_spikes[bins*i + j] == true) {
				add_int += pow(2, j);
			}
		}
		spikeInt.push_back(add_int);
	}
}

double BoolHistogram(vector<bool>& data) {
	int count = 0;
	for (vector<bool>::iterator it = data.begin(); it != data.end(); it++) {
		if (*it) count += 1;
	}
	return count * 1.0 / data.size();
}

void IntHistogram(vector<int>& data, vector<double> & histogram, int min, int max) {
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

void DoubleHistogram(vector<double> & data, vector<double> & histogram, double min, double max, double bin_width) {
	int number;
	number = ceil((max - min) / bin_width);
	vector<int> count(number, 0);
	int IND;
	for (vector<double>::iterator it = data.begin(); it != data.end(); it++) {
		IND = floor((*it - min) / bin_width);
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
		px = BoolHistogram(x);
		py = BoolHistogram(y);
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

		// Mutual information;
		double mi = 0;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				if (count_xy[i][j] != 0) {
					mi += count_xy[i][j] * 1.0 / np*log2(count_xy[i][j] * 1.0 / np / Px[i] / Py[j]);
				}
				else continue;
			}
		}
		if (mi<0) {
			cout << mi << endl;
			cout << "Press ENTER to continue:" << endl;
			cin.get();
		}
		return mi;
	}
}

double MI(vector<double>& x, vector<double>& y, double* x_max_and_min, double* y_max_and_min, double x_bin_size, double y_bin_size) {
	// total time period of X and Y;
	if (x.size() != y.size()) {
		cout << "Error: x and y don't have the same length." << endl;
		return 0;
	} else {
		// Calculate single probabilities;
		vector<double> Px, Py;
		DoubleHistogram(x, Px, x_max_and_min[1], x_max_and_min[0], x_bin_size);
		DoubleHistogram(y, Py, y_max_and_min[1], y_max_and_min[0], y_bin_size);

		// 	Calculate conditional probability;
		int time_bin_number = x.size(); // number of bins of time series;
		vector<int> add(Py.size(), 0);
		vector<vector<int> > count_xy(Px.size(), add);
		int indx, indy;
		for (int i = 0; i < time_bin_number; i++) {
			indx = floor((x[i] - x_max_and_min[1]) / x_bin_size);
			indy = floor((y[i] - y_max_and_min[1]) / y_bin_size);
			if (indx == Px.size()) indx --;
			if (indy == Py.size()) indy --;
			count_xy[indx][indy]++;
		}

		// Calculate mutual information;
		double mi = 0;
		double Pxy;
		for (int i = 0; i < Px.size(); i++) {
			for (int j = 0; j < Py.size(); j++) {
				if (count_xy[i][j] != 0) {
					Pxy = count_xy[i][j] * 1.0 / time_bin_number;
					mi += Pxy * log2(Pxy / (Px[i] * Py[j]));
				}
			}
		}
		if (mi<0) {
			cout << mi << endl;
			cout << "Press ENTER to continue:" << endl;
			cin.get();
		}
		return mi;
	}
}

double MI(vector<double>& x, vector<double>& y) {
	//	Compare the length of x & y;
	int length;
	if (x.size() > y.size()) length = y.size();
	else length = x.size();
	// Calculate actual occupancy;
	int bin_number = floor(sqrt(length / 5)); // bin_number: numberber of non-uniform bins;
	int occupancy = floor(length / bin_number);
	int x_residue = x.size() - bin_number*occupancy; // residue: remaining number of data that not used.
	int y_residue = y.size() - bin_number*occupancy;
	int data_size = bin_number*occupancy;

	// Find edges of histogram;
	vector<double> x_edges, y_edges;
	FindEdges(x, x_edges, occupancy, x_residue);
	FindEdges(y, y_edges, occupancy, y_residue);

	// 	Calculate conditional probability;
	vector<vector<double> > jointpdf;
	JointPDF(x, y, x_edges, y_edges, jointpdf);

	// Calculate mutual information;
	double mi = 0;
	for (int i = 0; i < bin_number; i++) {
		for (int j = 0; j < bin_number; j++) {
			if (jointpdf[i][j] != 0) {
				mi += jointpdf[i][j] * log2(pow(bin_number, 2)*jointpdf[i][j]);
			}
		}
	}
	if (mi<0) {
		cout << mi << endl;
		cout << "Press ENTER to continue:" << endl;
		cin.get();
	}
	return mi;
}

double MI(vector<double>& x, vector<double>& y, string pdf_path) {
	//	Compare the length of x & y;
	int length;
	if (x.size() > y.size()) length = y.size();
	else length = x.size();
	// Calculate actual occupancy;
	int bin_number = floor(sqrt(length / 5)); // bin_number: numberber of non-uniform bins;
	int occupancy = floor(length / bin_number);
	int x_residue = x.size() - bin_number*occupancy; // residue: remaining number of data that not used.
	int y_residue = y.size() - bin_number*occupancy;
	int data_size = bin_number*occupancy;

	// Find edges of histogram;
	vector<double> x_edges, y_edges;
	FindEdges(x, x_edges, occupancy, x_residue);
	FindEdges(y, y_edges, occupancy, y_residue);

	// 	Calculate conditional probability;
	vector<vector<double> > jointpdf;
	JointPDF(x, y, x_edges, y_edges, jointpdf);
	Print2D(pdf_path, "trunc", jointpdf);

	// Calculate mutual information;
	double mi = 0;
	for (int i = 0; i < bin_number; i++) {
		for (int j = 0; j < bin_number; j++) {
			if (jointpdf[i][j] != 0) {
				mi += jointpdf[i][j] * log2(pow(bin_number, 2)*jointpdf[i][j]);
			}
		}
	}
	if (mi<0) {
		cout << mi << endl;
		cout << "Press ENTER to continue:" << endl;
		cin.get();
	}
	return mi;
}

double MI(vector<bool>& binary_spikes, vector<double>& LFP, int time_bin_number, double LFP_max, double LFP_min, double bin_width) {
	// Prepare historgram;
	double p_spike;
	p_spike = BoolHistogram(binary_spikes);
	vector<double> LFP_histogram;
	DoubleHistogram(LFP, LFP_histogram, LFP_min, LFP_max, bin_width);

	// Calculate conditional probability;
	int histogram_bin_number = LFP_histogram.size();
	vector<int> count_xy_true;
	vector<int> count_xy_false;
	count_xy_true.resize(histogram_bin_number, 0);
	count_xy_false.resize(histogram_bin_number, 0);
	int ind;
	for (int i = 0; i < time_bin_number; i++) {
		ind = floor((LFP[i] - LFP_min) / bin_width);
		if (binary_spikes[i] == true) count_xy_true[ind]++;
		else count_xy_false[ind]++;
	}

	// Calculate mutual information;

	double mi = 0;
	for (int i = 0; i < histogram_bin_number; i++) {
		if (LFP_histogram[i] != 0) {
			if (count_xy_true[i] != 0) {
				mi += count_xy_true[i] * 1.0 / time_bin_number*log2(count_xy_true[i] * 1.0 / time_bin_number / p_spike / LFP_histogram[i]);
			}
			if (count_xy_false[i] != 0) {
				mi += count_xy_false[i] * 1.0 / time_bin_number*log2(count_xy_false[i] * 1.0 / time_bin_number / (1 - p_spike) / LFP_histogram[i]);
			}
		}
	}
	if (mi > 10000 or mi < 0) {
		cout << mi << endl;
		cout << "Press ENTER to continue:" << endl;
		cin.get();
	}
	return mi;
}

double MI(vector<bool>& binary_spikes, vector<double>& LFP, int bin_number) {
	// calculate the occupancy;
	// int bin_number = floor(sqrt(LFP.size() / expected_occupancy)); // bin_number: numberber of non-uniform bins;
	int occupancy = LFP.size() / bin_number;
	int residue = LFP.size() % bin_number; // residue: remaining number of data that not used.
	int data_size = LFP.size() - residue;

	// calculate histogram;
	vector<double> edges;
	FindEdges(LFP, edges, occupancy, residue);
	double p_spike;
	p_spike = BoolHistogram(binary_spikes);
	double pr_spikes[2]; // pr_spikes[0] = probability of spikes; pr_spikes[1] = probability of non-spikes;
	pr_spikes[0] = p_spike;
	pr_spikes[1] = 1 - p_spike;

	// Calculate conditional probability;
	vector<vector<double> > jointpdf;
	JointPDF(binary_spikes, LFP, edges, jointpdf);

	// Calculate mutual information;

	double mi = 0; // Pxy = P(x,y);
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j< bin_number; j++) {
			if (jointpdf[i][j] != 0) {
				mi += jointpdf[i][j] *log2(bin_number*jointpdf[i][j] / pr_spikes[i]);
			}
		}
	}
	if (mi<0) {
		cout << mi << endl;
		cout << "Press ENTER to continue:" << endl;
		cin.get();
	}
	return mi;
}

void TDMI(vector<double>& x, vector<double>& y, double dt, double tmax, int negative_time_delay, int positive_time_delay, vector<double> & tdmi) {
	vector<bool> x_bool, y_bool;
	ConvertSpikeToBinary(x, x_bool, tmax, dt);
	ConvertSpikeToBinary(y, y_bool, tmax, dt);
	tdmi.resize(positive_time_delay + negative_time_delay + 1, 0);

	// No shift;
	tdmi[negative_time_delay] = MI(x_bool, y_bool);

	// Negative shift;
	vector<bool> x_copy = x_bool, y_copy = y_bool;
	for (int i = 0; i < negative_time_delay; i++) {
		x_copy.erase(x_copy.begin());
		y_copy.erase(y_copy.end() - 1);
		tdmi[negative_time_delay - i - 1] = MI(x_copy, y_copy);
	}

	// Positive shift;
	x_copy = x_bool;
	y_copy = y_bool;
	for (int i = 0; i < positive_time_delay; i++) {
		x_copy.erase(x_copy.end() - 1);
		y_copy.erase(y_copy.begin());
		tdmi[negative_time_delay + i + 1] = MI(x_copy, y_copy);
	}
}

void TDMI_uniform(vector<double>& x, vector<double>& y, int negative_time_delay, int positive_time_delay, vector<double> & tdmi) {
	int repeat_number = positive_time_delay + negative_time_delay + 1;
	tdmi.resize(repeat_number, 0);
	// Measure the boundaries;
	int length = x.size();
	int numofpart = floor(sqrt(length));
	double x_max_and_min[2], y_max_and_min[2];
	x_max_and_min[0] = Max(x);
	x_max_and_min[1] = Min(x);
	y_max_and_min[0] = Max(y);
	y_max_and_min[1] = Min(y);
	double x_bin_size = (x_max_and_min[0] - x_max_and_min[1]) / numofpart;
	double y_bin_size = (y_max_and_min[0] - y_max_and_min[1]) / numofpart;

	// No shift;
	tdmi[negative_time_delay] = MI(x, y, x_max_and_min, y_max_and_min, x_bin_size, y_bin_size);

	// Negative shift;
	vector<double> x_copy = x, y_copy = y;
	for (int i = 0; i < negative_time_delay; i++) {
		x_copy.erase(x_copy.begin(), x_copy.begin() + 1);
		y_copy.erase(y_copy.end() - 1, y_copy.end());
		tdmi[negative_time_delay - i - 1] = MI(x_copy, y_copy, x_max_and_min, y_max_and_min, x_bin_size, y_bin_size);
	}

	// Positive shift;
	x_copy = x;
	y_copy = y;
	for (int i = 0; i < positive_time_delay; i++) {
		x_copy.erase(x_copy.end() - 1, x_copy.end());
		y_copy.erase(y_copy.begin(), y_copy.begin() + 1);
		tdmi[negative_time_delay + i + 1] = MI(x_copy, y_copy, x_max_and_min, y_max_and_min, x_bin_size, y_bin_size);
	}
}

void TDMI_adaptive(vector<double>& x, vector<double>& y, int negative_time_delay, int positive_time_delay, vector<double> & tdmi) {
	int repeat_number = positive_time_delay + negative_time_delay + 1;
	tdmi.resize(repeat_number, 0);

	string pdf_path;
	int index = 0;
	// Negative shift;
	vector<double> x_copy = x, y_copy = y;
	for (int i = 0; i < negative_time_delay; i++) {
		x_copy.erase(x_copy.begin(), x_copy.begin() + 1);
		y_copy.erase(y_copy.end() - 1, y_copy.end());
		pdf_path.clear();
		int2str(index, pdf_path);
		pdf_path = "../../results/jointpdfs/jpdf_" + pdf_path;
		index ++;
		tdmi[negative_time_delay - i - 1] = MI(x_copy, y_copy, pdf_path);
	}

	// No shift;
	pdf_path.clear();
	int2str(index, pdf_path);
	pdf_path = "../../results/jointpdfs/jpdf_" + pdf_path;
	index ++;
	tdmi[negative_time_delay] = MI(x, y, pdf_path);

	// Positive shift;
	x_copy = x;
	y_copy = y;
	for (int i = 0; i < positive_time_delay; i++) {
		x_copy.erase(x_copy.end() - 1, x_copy.end());
		y_copy.erase(y_copy.begin(), y_copy.begin() + 1);
		pdf_path.clear();
		int2str(index, pdf_path);
		pdf_path = "../../results/jointpdfs/jpdf_" + pdf_path;
		index ++;
		tdmi[negative_time_delay + i + 1] = MI(x_copy, y_copy, pdf_path);
	}
}

void TDMI(vector<double>& spikes, vector<double>& LFP, double dt, double sampling_dt, double bin_width, int negative_time_delay, int positive_time_delay, vector<double>& tdmi) {
	int repeat_number = negative_time_delay + positive_time_delay + 1;
	tdmi.resize(repeat_number, 0);
	double tmax = LFP.size()*sampling_dt;
	vector<bool> binary_spikes;
	ConvertSpikeToBinary(spikes, binary_spikes, tmax, dt);

	//random_shuffle(binary_spikes.begin(), binary_spikes.end());

	//	Calculate the mean local field potential;
	int time_bin_number = binary_spikes.size();
	int n = dt / sampling_dt; // number of LFP data point in single time step;
	vector<double> mean_LFP(time_bin_number, 0);
	for (int i = 0; i < time_bin_number; i++) {
		for (int j = 0; j < n; j++) mean_LFP[i] += LFP[i*n + j];
		mean_LFP[i] /= n;
	}

	// Measure Range of data;
	double LFP_max_and_min[2];
	FindMaxMin(mean_LFP, LFP_max_and_min, bin_width);

	// No shift;
	tdmi[negative_time_delay] = MI(binary_spikes, mean_LFP, time_bin_number, LFP_max_and_min[0], LFP_max_and_min[1], bin_width);

	// Negative shift;
	vector<bool> spike_copy = binary_spikes;
	vector<double> LFP_copy = mean_LFP;
	for (int i = 0; i < negative_time_delay; i++) {
		spike_copy.erase(spike_copy.begin());
		LFP_copy.erase(LFP_copy.end() - 1, LFP_copy.end());
		tdmi[negative_time_delay - i - 1] = MI(spike_copy, LFP_copy, time_bin_number - 1, LFP_max_and_min[0], LFP_max_and_min[1], bin_width);
	}

	// Positive shift;
	spike_copy = binary_spikes;
	LFP_copy = mean_LFP;
	for (int i = 0; i < positive_time_delay; i++) {
		spike_copy.erase(spike_copy.end() - 1);
		LFP_copy.erase(LFP_copy.begin(), LFP_copy.begin() + 1);
		tdmi[negative_time_delay + i + 1] = MI(spike_copy, LFP_copy, time_bin_number - 1, LFP_max_and_min[0], LFP_max_and_min[1], bin_width);
	}
}


void TDMI(vector<double>& spikes, vector<double>& LFP, double dt, double sampling_dt, int negative_time_delay, int positive_time_delay, vector<double>& tdmi, bool random_switch) {
	int dn = dt / sampling_dt; // number of LFP data points within single time step;
	int time_bin_number = floor(LFP.size() / dn); // number of reduced LFP data point;

	int bin_number = floor(spikes.size() / 5); // expected_occupancy = 5;
	// int bin_number = floor(time_bin_number / expected_occupancy / 2);

	// rearrange LFP data;
	vector<double> mean_LFP(time_bin_number, 0);
	double mean;
	for (int i = 0; i < time_bin_number; i++) {
		for (int j = 0; j < dn; j++) mean_LFP[i] += LFP[i*dn + j];
		mean_LFP[i] /= dn;
	}

	// prepare the spiking train;
	double tmax = mean_LFP.size() * dt;
	vector<bool> binary_spikes;
	ConvertSpikeToBinary(spikes, binary_spikes, tmax, dt);

	if (random_switch == true) {
		random_shuffle(binary_spikes.begin(), binary_spikes.end());
	}

	int repeat_number = negative_time_delay + positive_time_delay + 1;
	tdmi.resize(repeat_number, 0);
	int progress_counter = 0;
	double progress;
	char cr = (char)13;
	// No shift;
	tdmi[negative_time_delay] = MI(binary_spikes, mean_LFP, bin_number);
	progress_counter ++;
	progress = progress_counter*100.0/repeat_number;
	cout << cr << ">> Processing ... ";
	printf("%6.2f",progress);
	cout << "%";


	// Negative shift;
	vector<bool> spikes_copy = binary_spikes;
	vector<double> LFP_copy = mean_LFP;
	for (int i = 0; i < negative_time_delay; i++) {
		spikes_copy.erase(spikes_copy.begin());
		LFP_copy.erase(LFP_copy.end() - 1);
		tdmi[negative_time_delay - i - 1] = MI(spikes_copy, LFP_copy, bin_number);
		progress_counter ++;
		progress = progress_counter*100.0/repeat_number;
		cout << cr << ">> Processing ... ";
		printf("%6.2f",progress);
		cout << "%";
	}

	// Positive shift;
	spikes_copy = binary_spikes;
	LFP_copy = mean_LFP;
	for (int i = 0; i < positive_time_delay; i++) {
		spikes_copy.erase(spikes_copy.end() - 1);
		LFP_copy.erase(LFP_copy.begin());
		tdmi[negative_time_delay + i + 1] = MI(spikes_copy, LFP_copy, bin_number);
		progress_counter ++;
		progress = progress_counter*100.0/repeat_number;
		cout << cr << ">> Processing ... ";
		printf("%6.2f",progress);
		cout << "%";
	}
	cout << endl;
}
