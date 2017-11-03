//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-06-03
//	Description: Mutual information analysis program; version 1.0
//***************
#include "../include/mi_uniform.h"
#include "../include/io.h"
#include "../include/vecmanip.h"
#include <iostream>
#include <fstream>
#include <sstream>
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
	vector<bool> bool_series;
	vector<double> double_series;
	Read1D(argv[1], bool_series, 0, 1);
	Read1D(argv[2], double_series, 0, 1);
	double autoscale = 20; // with unit millisecond;
	double dt = 0.5;
	size_t length = floor(autoscale / dt);
	cout << length << '\n';
	size_t N = floor(bool_series.size() / length);
	if (bool_series.begin() + length*N != bool_series.end()) {
		bool_series.erase(bool_series.begin() + length*N, bool_series.end());
		double_series.erase(double_series.begin() + length*N, double_series.end());
	}
	size_t shape[2];
	shape[0] = N;
	shape[1] = length;
	vector<vector<bool> > boolmat, newboolmat;
	vector<vector<double> > doublemat, newdoublemat;
	Reshape(bool_series, boolmat, shape);
	Reshape(double_series, doublemat, shape);
	Transpose(boolmat, newboolmat);
	Transpose(doublemat, newdoublemat);
	// Set time range;
	size_t indx = atoi(argv[3]);
	istringstream range_in(argv[4]);
	string buffer;
	getline(range_in, buffer, ',');
	int negative_time_delay = atoi(buffer.c_str());
	getline(range_in, buffer, ',');
	int positive_time_delay = atoi(buffer.c_str());

	vector<bool> bool_copy = newboolmat[indx];
	vector<vector<double> > double_copy(newdoublemat.begin() + indx - negative_time_delay, newdoublemat.begin() + indx + positive_time_delay + 1);

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
	cout << "It takes " << (finish - start)*1.0 / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
