//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-06-03
//	Description: Mutual information analysis program; version 1.0
//***************
#include "../include/mi_uniform.h"
#include "../include/io.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using namespace std;
// compact function to calculate mutual information between multi-type signal
//	arguments:
//	argv[1] = path for bool series;
//	argv[2] = path for double series;
//  argv[3] = index of x variable;
//	argv[4] = range of timelag;
//	argv[5] = bin number of pdf of the double variable;
int main(int argc, const char* argv[]) {
	if (argc != 6) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// INPUT NEURONAL DATA:
	vector<vector<bool> > bool_series;
	vector<vector<double> > double_series;
	Read2D(argv[1], bool_series);
	Read2D(argv[2], double_series);
	// Set time range;
	size_t indx = atoi(argv[3]);
	string range = argv[4];
	string::size_type pos = range.find_first_of(',', 0 );
	int negative_time_delay = atoi(range.substr(0, pos).c_str());
	range.erase(0, pos + 1);
	int positive_time_delay = atoi(range.c_str());
	range = "";

	vector<bool> bool_copy = bool_series[indx];
	vector<vector<double> > double_copy(double_series.begin() + indx - negative_time_delay, double_series.begin() + indx + positive_time_delay + 1);

	size_t bin_num = atoi(argv[5]);
	vector<double> tdmi;
	TDMI(bool_copy, double_copy, tdmi, bin_num);

	//	Output data:
	ofstream data_out;
	cout << ">> Outputing data ... " << endl;
	data_out.open("./data/mi/mi_bd.csv");
	data_out << "timelag,mi" << endl;
	for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i++) {
		data_out << (i - negative_time_delay) << ',' << (double)tdmi[i] << '\n';
	}
	data_out.close();

	finish = clock();
	// Time counting:
	double ToTtime;
	ToTtime = (finish - start) / CLOCKS_PER_SEC;
	cout << "It takes " << (double)ToTtime << "s" << endl;
	return 0;
}
