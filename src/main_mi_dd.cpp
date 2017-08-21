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
//	argv[1] = mode code;
//	argv[2] = timing step for TDMI;
//	argv[3] = delay time range;
int main(int argc, const char* argv[]) {
	if (argc != 4) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// Preparing input args;
	vector<double> double_series_1, double_series_2;
	Read1D(argv[1], 0, 1, double_series_1);
	Read1D(argv[2], 0, 1, double_series_2);
  double dt = atof(argv[2]);
  string range = argv[3];
  string::size_type pos = range.find_first_of(',', 0 );
  int negative_time_delay = atoi(range.substr(0, pos).c_str());
  range.erase(0, pos + 1);
  int positive_time_delay = atoi(range.c_str());
  range = "";

  vector<double> tdmi;
	cout << ">> Calculating ordered TDMI ... " << endl;
	TDMI_adaptive(double_series_1, double_series_2, negative_time_delay, positive_time_delay, tdmi);

  // Output data:
	ofstream data_out;
	cout << ">> Outputing data ... " << endl;
	data_out.open("./data/mi/mi_dd.csv");
	data_out << "timelag,mi" << endl;
	for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i++) {
		data_out << (i - negative_time_delay) << ',' << (double)tdmi[i];
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
