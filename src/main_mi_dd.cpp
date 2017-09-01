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
#include <stdexcept>

using namespace std;
// compact function to calculate mutual information between multi-type signal
//	arguments:
//	argv[1] = path for time series x, m by n, indicate m variables with n trials;
//	argv[2] = path for time series y;
//	argv[3] = index of x variable in series;
//	argv[4] = delayed time range, seperated by comma;
int main(int argc, const char* argv[]) {
	if (argc != 5) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// Preparing input args;
	vector<vector<double> > double_series_1, double_series_2;
	Read2D(argv[1], double_series_1);
	Read2D(argv[2], double_series_2);
	int indx = atoi(argv[3]);
  string range = argv[4];
  int negative_time_delay = atoi(range.c_str());
	string::size_type pos = range.find_first_of(',', 0 );
  int positive_time_delay = atoi(range.c_str() + pos + 1);
  range = "";


	for (vector<vector<double> >::iterator it = double_series_1.begin(); it != double_series_1.end(); it ++) {
		for (vector<double>::iterator itt = it->begin(); itt != it->end(); itt ++) {
			cout << *itt << '\t';
		}
		cout << endl;
	}

	cout << endl;
	cout << endl;

	for (vector<vector<double> >::iterator it = double_series_2.begin(); it != double_series_2.end(); it ++) {
		for (vector<double>::iterator itt = it->begin(); itt != it->end(); itt ++) {
			cout << *itt << '\t';
		}
		cout << endl;
	}
	cout << endl;
	cout << endl;

	vector<double> s1 = double_series_1[indx];

	for (vector<double>::iterator it = s1.begin(); it != s1.end(); it ++) {
		cout << *it << '\t';
	}
	cout << endl;
	cout << s1.size() << endl;



	vector<vector<double> > s2(double_series_2.begin() + indx - negative_time_delay, double_series_2.begin() + indx + positive_time_delay + 1);
	for (vector<vector<double> >::iterator it = s2.begin(); it != s2.end(); it ++) {
		for (vector<double>::iterator itt = it->begin(); itt != it->end(); itt ++) {
			cout << *itt << '\t';
		}
		cout << endl;
	}
	cout << s2.size()	<< ',' << s2[0].size() << '\n';

	vector<double> tdmi;
	cout << ">> Calculating TDMI ... " << endl;
	TDMI(s1, s2, tdmi);

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
