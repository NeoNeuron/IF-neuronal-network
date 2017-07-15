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
#include <stdexcept>
using namespace std;

//	Function of calculating LFP with point current source model in 1-D loop network case;
//	arguments:
//	argv[1] = loading path of neural data;
//	argv[2] = list of indices of connected neurons, seperated by comma, that contribute to LFP;
//	argv[3] = time range, seperated by comma, with unit in milliseconds;
//	argv[4] = filename of output raster file;
int main(int argc, const char* argv[]) {
	if (argc != 5) {
		throw runtime_error("wrong number of args");
	}
	//	Defined folder path;
	string path = argv[1];
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
	printf(">> There are %d connected neuron which contribute to LFP.\n", neuron_num);

	//	Choose objective time range;
	double t_range[2]; // t_range[0] = t_min; t_range[1] = t_max;
	string range_str = argv[3];
  pos = range_str.find_first_of(',', 0);
  t_range[0] = atof(range_str.substr(0, pos).c_str());
  range_str.erase(0, pos + 1);
  t_range[1] = atof(range_str.c_str());
  range_str = "";
	printf(">> Time range is (%.2f, %.2f] ms.\n", t_range[0], t_range[1]);

	cout << ">> Calculating LFP ..." << endl;
	int total_neuron_number = atoi(argv[6]);
	vector<double> lfp;
	// LFP(t_range, total_neuron_number, list, potential_path, excitatory_conductance_path, inhibitory_conductance_path, lfp);
	LFP(t_range, list, path, lfp);

	//	Output lfp:
	cout << ">> Outputing LFP ..." << endl;
	string out_dir = "./data/lfp/";
	string lfp_path = out_dir + argv[4];
	OutLFP(lfp_path, lfp);
	cout << ">> Finished." << endl;
	return 0;
}
