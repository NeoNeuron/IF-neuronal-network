//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-13 15:07:34
//	Description: main file for lfp.h and lfp.cpp
//***************
#include "../include/lfp.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
using namespace std;

//	Function of calculating LFP with point current source model in 1-D loop network case;
//	arguments:
//	argv[1] = path of neural data;
//	argv[2] = list of indices of connected neurons, seperated by comma, that contribute to LFP;
//	argv[3] = time range, seperated by comma, with unit in milliseconds;
//  argv[4] = binning size of binary time series of spike train;
//	argv[5] = path of output raster file;
int main(int argc, const char* argv[]) {
	if (argc != 6) {
		throw runtime_error("wrong number of args");
	}
	cout << "==========" << endl;
	//	Defined folder path;
	//	Analyze listing series;
	string list_str =  argv[2];
	vector<int> list;
	string ss;
	string::size_type pos = list_str.find_first_of(',', 0);
	while (pos != list_str.npos) {
		ss = list_str.substr(0, pos);
		list.push_back(atoi(ss.c_str()));
		list_str.erase(0, pos + 1);
		ss.clear();
		pos = list_str.find_first_of(',', 0);
	}
	if (list_str.size() != 0) {
		list.push_back(atoi(list_str.c_str()));
	}
	sort(list.begin(),list.end());
	int neuron_num = list.size();
	printf(">> %d connected neuron contribute to LFP\n", neuron_num);

	//	Choose objective time range;
	double t_range[2]; // t_range[0] = t_min; t_range[1] = t_max;
	string range_str = argv[3];
  pos = range_str.find_first_of(',', 0);
  t_range[0] = atof(range_str.substr(0, pos).c_str());
  range_str.erase(0, pos + 1);
  t_range[1] = atof(range_str.c_str());
  range_str = "";
	printf(">> Time range = (%.2f, %.2f] ms\n", t_range[0], t_range[1]);

	cout << ">> Calculating LFP ..." << endl;
	vector<double> lfp;
	LFP(t_range, list, argv[1], lfp);

	// rearrange LFP data;
	double dt_sampling = 0.03125;
	int dn = atof(argv[4]) / dt_sampling; // number of LFP data points within single time step;
	int num_bin = floor((t_range[1] - t_range[0]) / atof(argv[4])); // number of reduced LFP data point;

	vector<double> mean_lfp(num_bin, 0);
	for (int i = 0; i < num_bin; i++) {
		for (int j = 0; j < dn; j++) mean_lfp[i] += lfp[i*dn + j];
		mean_lfp[i] /= dn;
	}

	//	Output lfp:
	cout << ">> Outputing LFP ..." << endl;
	OutLFP(argv[5], mean_lfp);
	cout << ">> Done" << endl << "==========" << endl;
	return 0;
}
