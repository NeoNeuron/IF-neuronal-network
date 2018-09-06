// ***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2018-09-02
//	Description: main file for lfp.h and lfp.cpp
//***************
#include "../include/lfp.h"
#include "../include/io.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
using namespace std;

//	Function of calculating LFP with point current source model in 1-D loop network case;
//
//	arguments:
//
//	argv[1] = path of neural data;
//	argv[2] = path of output LFP file;
//	argv[3] = list of indices of connected neurons, seperated by comma, that contribute to LFP;
//	argv[4] = time range, seperated by comma, with unit in milliseconds;
//	argv[5] = binning size of binary time series of spike train;
//
int main(int argc, const char* argv[]) {
	if (argc != 6) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// Analyze listing series;
	vector<int> list;
	stringstream inlist(argv[3]);
	string buffer;
	while (getline(inlist, buffer, ',')) {
		list.push_back(atoi(buffer.c_str()));
	}
	int neuron_num = list.size();
	printf(">> %d connected neuron contribute to LFP\n", neuron_num);

	//	Choose objective time range;
	double t_range[2]; // t_range[0] = t_min; t_range[1] = t_max;
  stringstream irange(argv[4]);
  for (size_t i = 0; i < 2; i++) {
    getline(irange, buffer, ',');
    t_range[i] = atof(buffer.c_str());
  }
	printf(">> Time range = (%.2f, %.2f] ms\n", t_range[0], t_range[1]);

	vector<double> lfp;
	LFP(argv[1], lfp, list, t_range);

	//	Output lfp:
	Print1DBin(argv[2], lfp, "trunc");
	finish = clock();
	// counting time;
	cout << "[-] LFP generation takes " << (finish - start)*1.0 / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
