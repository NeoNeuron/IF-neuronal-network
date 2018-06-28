//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2018-01-27
//	Description: Mutual information analysis program; version 1.0
//***************
#include "../include/mi_uniform.h"
#include "../include/io.h"
#include "../include/vecmanip.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <pthread.h>
#include <stdexcept>

using namespace std;
#define THREADS_NUM = 5;
// compact function to calculate mutual information between multi-type signal
//	arguments:
//	argv[1] = path for bool series;
//	argv[2] = path for double series;
//	argv[3] = path of output mi file;
//	argv[4] = range of timelag;
//	argv[5] = bin size of pdf of the double variable;
struct thread_args {
	vector<bool> bool_series;
	vector<double> double_series;
	double binsize;
};

double 
void TDMI(vector<bool>& x, vector<double>& y, vector<double>& tdmi, size_t* range, double binsize) {
	// initialize container of tdmi;
	tdmi.clear();
	tdmi.resize(range[0] + range[1] + 1, 0);
	// prepare series;
	size_t res = range[0];
	if (res < range[1]) res = range[1];
	vector<bool> x_copy(x.begin(), x.end() - res);
	vector<double> y_copy(y.begin(), y.end() - res);
	// No shift;
	tdmi[range[0]] = MIBD(x_copy, y_copy, binsize, false);
	// Negative shift;
	for (size_t i = 0; i < range[0]; i++) {
		x_copy.erase(x_copy.begin());
		x_copy.insert(x_copy.end(), *(x.end() - res + i));
		tdmi[range[0] - i - 1] = MIBD(x_copy, y_copy, binsize, false);
	}
	// Positive shift;
	x_copy.clear();
	x_copy.insert(x_copy.end(), x.begin(), x.end() - res);
	for (size_t i = 0; i < range[1]; i++) {
		y_copy.erase(y_copy.begin());
		y_copy.insert(y_copy.end(), *(y.end() - res + i));
		tdmi[range[0] + i + 1] = MIBD(x_copy, y_copy, binsize, false);
	}
	return;
}
int main(int argc, const char* argv[]) {
	if (argc != 6) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// INPUT NEURONAL DATA:
	vector<bool> bool_series;
	vector<double> double_series;
	string buffer;
	Read1DBin(argv[1], bool_series, 0, 0);
	Read1DBin(argv[2], double_series, 0, 0);
	// Set time range;
	size_t range[2];
	istringstream range_in(argv[4]);
	getline(range_in, buffer, ',');
	int ntd = atoi(buffer.c_str());
	range[0] = ntd;
	getline(range_in, buffer, ',');
	int ptd = atoi(buffer.c_str());
	range[1] = ptd;
	// Calculate mutual information;
	double binsize = atof(argv[5]);
	vector<double> tdmi;
	TDMI(bool_series, double_series, tdmi, range, binsize);

	//	Output data:
	ofstream data_out;
	cout << ">> Outputing data ... " << endl;
	data_out.open(argv[3]);
	data_out << "timelag,mi" << endl;
	for (int i = 0; i < ntd + ptd + 1; i++) {
		data_out << i - ntd << ',' << setprecision(15) << (double)tdmi[i] << '\n';
	}
	data_out.close();

	finish = clock();
	// Time counting:
	cout << "It takes " << (finish - start)*1.0 / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
