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
//	argv[5] = bin number of pdf of the double variable;
//	argv[6] = mutual info calculation algorithms; 'direct', 'partial' or 'full';
int main(int argc, const char* argv[]) {
	if (argc != 7) throw runtime_error("wrong number of args");
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
	size_t bin_num = atoi(argv[5]);
	vector<double> tdmi;
	if (strcmp(argv[6], "direct") == 0) {
		TDMI(bool_series, double_series, tdmi, range, bin_num);
	} else {
		// autocovariance operations;
		double autoscale = 20; // with unit millisecond;
		double dt = atof(argv[4]);
		size_t length = floor(autoscale / dt);
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
		// Print2D("bool_test", newboolmat, "trunc");
		// Print2D("double_test", newdoublemat, "trunc");
		if (strcmp(argv[6], "full") == 0) {
			TDMI(newboolmat, newdoublemat, tdmi, range, bin_num);
		} else if (strcmp(argv[6], "partial") == 0) {
			size_t indy = 20;
			vector<vector<bool> > bool_copy(newboolmat.begin() + indy - ptd, newboolmat.begin() + indy + ntd + 1);
			vector<double> double_copy = newdoublemat[indy];
			TDMI(bool_copy, double_copy, tdmi, bin_num);
		} else throw runtime_error("incorrect algorithm choosing");
	}

	//	Output data:
	ofstream data_out;
	cout << ">> Outputing data ... " << endl;
	data_out.open("./data/mi/mi_bd.csv");
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
