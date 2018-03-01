//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-06-03
//	Description: Mutual information analysis program; version 1.0
//***************
#include "../include/mi_uniform.h"
#include "../include/io.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <stdexcept>

using namespace std;
int myrandom(int i) {return rand()%i;}
// compact function to calculate mutual information between multi-type signal
//	arguments:
//	argv[1] = path for time series x, m by n, indicate m variables with n trials;
//	argv[2] = path for time series y;
//	argv[3] = index of x variable in series;
//	argv[4] = delayed time range, seperated by comma;
//  argv[5] = bin number of variable x;
//  argv[6] = bin number of variable y;
int main(int argc, const char* argv[]) {
	if (argc != 7) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// Preparing input args;
	vector<vector<double> > double_series_1, double_series_2;
	Read2D(argv[1], double_series_1);
	Read2D(argv[2], double_series_2);
	int indx = atoi(argv[3]);
	istringstream range_in(argv[4]);
	string buffer;
	getline(range_in, buffer, ',');
	int negative_time_delay = atoi(buffer.c_str());
	getline(range_in, buffer, ',');
	int positive_time_delay = atoi(buffer.c_str());

	vector<double> s1 = double_series_1[indx];
	vector<vector<double> > s2(double_series_2.begin() + indx - negative_time_delay, double_series_2.begin() + indx + positive_time_delay + 1);

	size_t x_bin_num = atoi(argv[5]), y_bin_num = atoi(argv[6]);
	vector<double> tdmi;
	TDMI(s1, s2, tdmi, x_bin_num, y_bin_num);
	// // Randomly shuffle double_series_2;
	// srand(unsigned(time(0)));
	// random_shuffle(double_series_2.begin(), double_series_2.end(), myrandom);
	// vector<vector<double> > s2_shuffle(double_series_2.begin() + indx - negative_time_delay, double_series_2.begin() + indx + positive_time_delay + 1);
	// vector<double> tdmi_shuffle;
	// TDMI(s1, s2_shuffle, tdmi_shuffle, x_bin_num, y_bin_num);


  // Output data:
	ofstream data_out;
	data_out.open("./data/mi/mi_dd.csv");
	data_out << "timelag,mi" << endl;
	for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i++) {
		data_out << (i - negative_time_delay) << ',' << tdmi[i] << endl;
	}
	data_out.close();

	finish = clock();
	// Time counting:
	cout << "It takes " << (finish - start) * 1.0 / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
