//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-06-03
//	Description: Mutual information analysis program; version 1.0
//***************
#include "../include/mi.h"
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
//	argv[3] = range of timelag;
int main(int argc, const char* argv[]) {
	if (argc != 4) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// INPUT NEURONAL DATA:
	vector<bool> bool_series;
	vector<double> double_series;
	Read1D(argv[1], 0, 1, bool_series);
	Read1D(argv[2], 0, 1, double_series);
	// Set time range;
	string range = argv[3];
	string::size_type pos = range.find_first_of(',', 0 );
	int negative_time_delay = atoi(range.substr(0, pos).c_str());
	range.erase(0, pos + 1);
	int positive_time_delay = atoi(range.c_str());
	range = "";
	// int number = 0;
	// for (vector<bool>::iterator it = bool_series.begin(); it != bool_series.end(); it ++) {
	// 	if (*it) number ++;
	// }
	// cout << number << endl;
	vector<double> tdmi_ordered;
	vector<double> tdmi_random;
	cout << ">> Calculating ordered TDMI ... " << endl;
	TDMI(bool_series, double_series, negative_time_delay, positive_time_delay, tdmi_ordered, false);
	cout << ">> Calculating swapped TDMI ... " << endl;
	TDMI(bool_series, double_series, negative_time_delay, positive_time_delay, tdmi_random, true);

	//	Output data:
	ofstream data_out;
	cout << ">> Outputing data ... " << endl;
	data_out.open("./data/mi/mi_bd.csv");
	data_out << "timelag,mi,shuffle" << endl;
	for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i++) {
		data_out << (i - negative_time_delay) << ',' << (double)tdmi_ordered[i] << ',' << (double)tdmi_random[i];
		if (i < positive_time_delay + negative_time_delay) data_out << endl;
	}
	data_out.close();

	finish = clock();
	// Time counting:
	double ToTtime;
	ToTtime = (finish - start) / CLOCKS_PER_SEC;
	cout << "It takes " << (double)ToTtime << "s" << endl;
	return 0;
}
