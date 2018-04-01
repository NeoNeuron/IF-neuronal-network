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
#include <ctime>
#include <stdexcept>
using namespace std;
// compact function to calculate mutual information between multi-type signal
//	arguments:
//	argv[1] = path for time series x;
//	argv[2] = path for time series y;
//	argv[3] = delayed time range, seperated by comma;
//  argv[4] = bin size of histogram;
int main(int argc, const char* argv[]) {
	if (argc != 5) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// Preparing input args;
	vector<double> double_series_1, double_series_2;
	Read1DBin(argv[1], double_series_1, 0, 0);
	Read1DBin(argv[2], double_series_2, 0, 0);
	istringstream range_in(argv[3]);
	size_t range[2];
	string buffer;
	getline(range_in, buffer, ',');
	int negative_time_delay = atoi(buffer.c_str());
	range[0] = negative_time_delay;
	getline(range_in, buffer, ',');
	int positive_time_delay = atoi(buffer.c_str());
	range[1] = positive_time_delay;

	double binsize = atof(argv[4]);
	vector<double> tdmi;
	TDMI(double_series_1, double_series_2, tdmi, range, binsize);

  // Output data:
	ofstream data_out;
	data_out.open("./data/mi/mi_dd_LFP.csv");
	data_out << "timelag,mi" << endl;
	for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i++) {
		data_out << (i - negative_time_delay) << ',' << setprecision(15) << tdmi[i] << endl;
	}
	data_out.close();

	finish = clock();
	// Time counting:
	cout << "It takes " << (finish - start) * 1.0 / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
