//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-11-08
//	Description: Mutual information analysis program, between spikes and spikes;
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
//	argv[1] = path for first spike train;
//	argv[2] = path for second spike train;
//	argv[3] = range of timelag;
//	argv[4] = size of timing step;
int main(int argc, const char* argv[]) {
	if (argc != 5) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// INPUT NEURONAL DATA:
	vector<bool> x, y;
	Read1D(argv[1], x, 0, 1);
	Read1D(argv[2], y, 0, 1);
	// Set time range;
	size_t range[2];
	istringstream range_in(argv[3]);
	string buffer;
	getline(range_in, buffer, ',');
	int ntd = atoi(buffer.c_str());
	getline(range_in, buffer, ',');
	int ptd = atoi(buffer.c_str());
  range[0] = ntd;
  range[1] = ptd;
  double dt = atof(argv[4]);

	vector<double> tdmi;
	TDMI(x, y, tdmi, range);

	//	Output data:
	ofstream data_out;
	cout << ">> Outputing data ... " << endl;
	data_out.open("./data/mi/mi_bb.csv");
	data_out << "timelag,mi" << endl;
	for (int i = 0; i < ntd + ptd + 1; i++) {
		data_out << i - ntd << ',' << (double)tdmi[i] << '\n';
	}
	data_out.close();

	finish = clock();
	// Time counting:
	cout << "It takes " << (finish - start)*1.0 / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
