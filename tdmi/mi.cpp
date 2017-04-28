//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-02-27 21:48:11
//	Description: source file of mi.h
//***************
#include "mi.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

const double PI=3.1415926;

bool comp(const double &a, const double &b) {
	return a < b;
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
	edges.clear();
	for (int i = 0; i < N; i++) {
		edges.push_back(data_copy[i*occupancy]);
	}
	data_copy.clear();
}

void FindMaxMin(vector<double>& data, double *max_and_min, double bin_width) {
	double min, max;
	for (vector<double>::iterator it = data.begin(); it != data.end(); it++) {
		if (it == data.begin()) {
			min = *it;
			max = *it;
		} else {
			if (*it < min) min = *it;
			else if (*it > max) max = *it;
			else continue;
		}
	}
	min = floor(min / bin_width)*bin_width;
	max = ceil(max / bin_width)*bin_width;
	*max_and_min = max;
	*(max_and_min + 1) = min;
}

void ReadData(string filename, vector<vector<double> >& data) {
	data.clear();
	ifstream ifile;
	ifile.open(filename.c_str());
	string s;
	vector<double> add_double;
	while (getline(ifile, s)) {
		add_double.clear();
		string::size_type pos = s.find_first_of('\t', 0);
		string ss;
		while (pos != s.npos) {
			ss = s.substr(0, pos);
			add_double.push_back(atof(ss.c_str()));
			s.erase(0, pos + 1);
			ss.clear();
			pos = s.find_first_of('\t', 0);
		}
		pos = s.find_first_of('\n', 0);
		if (pos == 0) continue;
		else {
			ss = s.substr(0, pos);
			add_double.push_back(atof(ss.c_str()));
		}
		data.push_back(add_double);
		s.clear();
	}
	ifile.close();
}

void ReadData(string filename, vector<vector<int> >& data) {
	data.clear();
	ifstream ifile;
	ifile.open(filename.c_str());
	string s;
	vector<int> add_int;
	while (getline(ifile, s)) {
		add_int.clear();
		string::size_type pos = s.find_first_of('\t', 0);
		while (pos != s.npos) {
			string ss = s.substr(0, pos);
			add_int.push_back(atof(ss.c_str()));
			s.erase(0, pos + 1);
			ss.clear();
			pos = s.find_first_of('\t', 0);
		}
		data.push_back(add_int);
		s.clear();
	}
	ifile.close();
}

void ReadData(string filename, vector<double> & data) {
	data.clear();
	ifstream ifile;
	ifile.open(filename.c_str());
	string s;
	double add_double;
	string::size_type pos;
	string ss;
	while (getline(ifile, s)) {
		pos = s.find_first_of('\n', 0);
		ss = s.substr(0, pos);
		add_double = atof(ss.c_str());
		data.push_back(add_double);
	}
	ifile.close();
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
		if (*it == true) count += 1;
	}
	double length = data.size();
	return count / length;
}

void IntHistogram(vector<int>& data, vector<double> & histogram, int min, int max) {
	int number = max - min + 1;
	vector<int> count(number, 0);
	int IND;
	for (vector<int>::iterator it = data.begin(); it != data.end(); it++) {
		IND = *it - min;
		count[IND] ++;
	}
	for (int i = 0; i < number; i++) {
		histogram.push_back(count[i] * 1.0 / data.size());
	}
}

void DoubleHistogram(vector<double> & data, vector<double> & histogram, double min, double max, double bin_width) {
	int number;
	number = ceil((max - min) / bin_width);
	vector<int> count(number, 0);
	int IND;
	for (vector<double>::iterator it = data.begin(); it != data.end(); it++) {
		IND = floor((*it - min) / bin_width);
		count[IND] ++;
	}
	histogram.clear();
	for (int i = 0; i < number; i++) {
		histogram.push_back(count[i] * 1.0 / data.size());
	}
}

double MI(vector<bool>& x, vector<bool>& y, int bin_size) {
	if (x.size() != y.size()) {
		cout << "ERROR: x and y don't have the same length." << endl;
		return 0;
	} else {
		// Convert bool vector to int vector;

		vector<int> x_int, y_int;
		ConvertBinaryToInt(x, x_int, bin_size);
		ConvertBinaryToInt(y, y_int, bin_size);

		// Histograms;
		int max_value = pow(2, bin_size);
		vector<double> Px, Py;
		IntHistogram(x_int, Px, 0, max_value - 1);
		IntHistogram(y_int, Py, 0, max_value - 1);
		
		vector<int> add(max_value, 0);
		vector<vector<int> > count_xy(max_value, add);
		for (int i = 0; i < x_int.size(); i++) {
			count_xy[x_int[i]][y_int[i]]++;
		}

		// Mutual information;
		double mi = 0;
		for (int i = 0; i < max_value; i++) {
			for (int j = 0; j < max_value; j++) {
				if (count_xy[i][j] != 0) {
					mi += count_xy[i][j] * 1.0 / x_int.size()*log2(count_xy[i][j] * 1.0 / x_int.size() / Px[i] / Py[j]);
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

double MI(vector<double>& x, vector<double>& y, double expected_occupancy) {
	//	Compare the length of x & y;
	int length;
	if (x.size() > y.size()) length = y.size();
	else length = x.size();
	// Calculate actual occupancy;
	int bin_number = floor(sqrt(length / expected_occupancy)); // bin_number: numberber of non-uniform bins;
	int occupancy = floor(length / bin_number);
	int x_residue = x.size() - bin_number*occupancy; // residue: remaining number of data that not used.
	int y_residue = y.size() - bin_number*occupancy;
	int data_size = bin_number*occupancy;

	// Find edges of histogram;
	vector<double> x_edges, y_edges;
	FindEdges(x, x_edges, occupancy, x_residue);
	FindEdges(y, y_edges, occupancy, y_residue);

	// 	Calculate conditional probability;
	vector<int> add(bin_number, 0);
	vector<vector<int> > count_xy(bin_number, add);
	int indx, indy;
	for (int i = 0; i < data_size; i++) {
		indx = -1, indy = -1;
		for (int j = 0; j < bin_number; j++) {
			if (x[i] < x_edges[j]) {indx = j - 1; break;}
			if (indx == -1) indx = bin_number - 1;
		}
		for (int j = 0; j < bin_number; j++) {
			if (y[i] < y_edges[j]) {indy = j - 1; break;}
			if (indy == -1) indy = bin_number - 1;
		}
		count_xy[indx][indy]++;
	}

	// Calculate mutual information;
	double mi = 0, Pxy;
	for (int i = 0; i < bin_number; i++) {
		for (int j = 0; j < bin_number; j++) {
			if (count_xy[i][j] != 0) {
				Pxy = count_xy[i][j] * 1.0 / data_size;
				mi += Pxy * log2(pow(bin_number, 2)*Pxy);
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

double MI(vector<bool>& binary_spikes, vector<double>& LFP, int expected_occupancy) {
	// calculate the occupancy;	
	int bin_number = floor(sqrt(LFP.size() / expected_occupancy)); // bin_number: numberber of non-uniform bins;
	int occupancy = LFP.size() / bin_number, residue = LFP.size() % bin_number; // residue: remaining number of data that not used.
	int data_size = LFP.size() - residue;

	// calculate histogram;
	vector<double> edges;
	FindEdges(LFP, edges, occupancy, residue);
	double p_spike;
	p_spike = BoolHistogram(binary_spikes);

	// Calculate conditional probability;
	vector<int> count_xy_true(bin_number, 0);
	vector<int> count_xy_false(bin_number, 0);
	
	int ind;
	for (int i = 0; i < data_size; i++) {
		ind = -1;
		for (int j = 0; j < bin_number; j++) {
			if (LFP[i] < edges[j]) {ind = j - 1; break;}
			if (ind == -1) ind = bin_number - 1;
		}
		if (binary_spikes[i] == true) count_xy_true[ind]++;
		else count_xy_false[ind]++;
	}

	// Calculate mutual information;

	double mi = 0, Pxy = 0; // Pxy = P(x,y);
	for (int i = 0; i < bin_number; i++) {	
		if (count_xy_true[i] != 0) {
			Pxy = count_xy_true[i] * 1.0 / data_size;
			mi += Pxy *log2(bin_number*Pxy / p_spike);
		}
		if (count_xy_false[i] != 0) {
			Pxy = count_xy_false[i] * 1.0 / data_size;
			mi += Pxy *log2(bin_number*Pxy / (1 - p_spike));
		}
	}
	if (mi<0) {
		cout << mi << endl;
		cout << "Press ENTER to continue:" << endl;
		cin.get();
	}
	return mi;
}

void TDMI(vector<double>& x, vector<double>& y, double dt, double tmax, int bin_size, int negative_time_delay, int positive_time_delay, vector<double> & tdmi) {
	vector<bool> x_bool, y_bool;
	ConvertSpikeToBinary(x, x_bool, tmax, dt);
	ConvertSpikeToBinary(y, y_bool, tmax, dt);
	int repeat_number = positive_time_delay + negative_time_delay + 1;
	tdmi.resize(repeat_number, 0);

	// No shift;
	tdmi[negative_time_delay] = MI(x_bool, y_bool, bin_size);

	// Negative shift;
	vector<bool> x_copy = x_bool, y_copy = y_bool;
	for (int i = 0; i < negative_time_delay; i++) {
		x_copy.erase(x_copy.begin());
		y_copy.erase(y_copy.end() - 1);
		tdmi[negative_time_delay - i - 1] = MI(x_copy, y_copy, bin_size);
	}

	// Positive shift;
	x_copy = x_bool;
	y_copy = y_bool;
	for (int i = 0; i < positive_time_delay; i++) {
		x_copy.erase(x_copy.end() - 1);
		y_copy.erase(y_copy.begin());
		tdmi[negative_time_delay + i + 1] = MI(x_copy, y_copy, bin_size);
	}
}

void TDMI(vector<double>& x, vector<double>& y, double x_bin_size, double y_bin_size, int negative_time_delay, int positive_time_delay, vector<double> & tdmi) {
	int repeat_number = positive_time_delay + negative_time_delay + 1;
	tdmi.resize(repeat_number, 0);
	// Measure the boundaries;

	double x_max_and_min[2], y_max_and_min[2];
	FindMaxMin(x, x_max_and_min, x_bin_size);
	FindMaxMin(y, y_max_and_min, y_bin_size);

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

void TDMI(vector<double>& x, vector<double>& y, int expected_occupancy, int negative_time_delay, int positive_time_delay, vector<double> & tdmi) {
	int repeat_number = positive_time_delay + negative_time_delay + 1;
	tdmi.resize(repeat_number, 0);

	// No shift;
	tdmi[negative_time_delay] = MI(x, y, expected_occupancy);
	
	// Negative shift;
	vector<double> x_copy = x, y_copy = y;
	for (int i = 0; i < negative_time_delay; i++) {
		x_copy.erase(x_copy.begin(), x_copy.begin() + 1);
		y_copy.erase(y_copy.end() - 1, y_copy.end());
		tdmi[negative_time_delay - i - 1] = MI(x_copy, y_copy, expected_occupancy);
	}

	// Positive shift;
	x_copy = x;
	y_copy = y;
	for (int i = 0; i < positive_time_delay; i++) {
		x_copy.erase(x_copy.end() - 1, x_copy.end());
		y_copy.erase(y_copy.begin(), y_copy.begin() + 1);		
		tdmi[negative_time_delay + i + 1] = MI(x_copy, y_copy, expected_occupancy);
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


void TDMI(vector<double>& spikes, vector<double>& LFP, int expected_occupancy, double dt, double sampling_dt, int negative_time_delay, int positive_time_delay, vector<double>& tdmi, bool random_switch) {
	int dn = dt / sampling_dt; // number of LFP data points within single time step;
	int time_bin_number = floor(LFP.size() / dn); // number of reduced LFP data point;
	
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
	tdmi[negative_time_delay] = MI(binary_spikes, mean_LFP, expected_occupancy);
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
		tdmi[negative_time_delay - i - 1] = MI(spikes_copy, LFP_copy, expected_occupancy);
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
		tdmi[negative_time_delay + i + 1] = MI(spikes_copy, LFP_copy, expected_occupancy);
		progress_counter ++;
		progress = progress_counter*100.0/repeat_number;
		cout << cr << ">> Processing ... ";
		printf("%6.2f",progress);
		cout << "%";
	}
	cout << endl;
}
