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
#include <stdexcept>

using namespace std;
// compact function to calculate mutual information between multi-type signal
//	arguments:
//	argv[1] = path for bool series;
//	argv[2] = path for double series;
//	argv[3] = range of timelag;
//	argv[4] = size of timing step;
//	argv[5] = threshold of pdf for the double variable;
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
	istringstream range_in(argv[3]);
	getline(range_in, buffer, ',');
	int ntd = atoi(buffer.c_str());
	range[0] = ntd;
	getline(range_in, buffer, ',');
	int ptd = atoi(buffer.c_str());
	range[1] = ptd;
	// Calculate mutual information;
	double threshold = atof(argv[5]);
	vector<double> tdmi;
	TDMI(bool_series, double_series, tdmi, range, threshold);

	//	Output data:
	ofstream data_out;
	cout << ">> Outputing data ... " << endl;
	data_out.open("./data/mi/mi_bd_2bins.csv");
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
