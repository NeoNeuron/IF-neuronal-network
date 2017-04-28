//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-13 15:07:34
//	Description: test file for lfp.h and lfp.cpp
//***************
#include "lfp.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>

using namespace std;

//	Function of calculating LFP with point current source model in 1-D loop network case;
//	arguments:
//	argv[1] = loading directory for neural data;
//	argv[2] = neuronal index for target neuron in pre network;
//	argv[3] = list of indices of connected neurons that contribute to LFP;
//	argv[4] = lower bond of time range;
//	argv[5] = upper bond of time range;
//	argv[6] = total neuron number;
int main(int argc, const char* argv[]) {
	if (argc != 7) {
		throw runtime_error("wrong number of args");
	}
	//	Defined folder path;
	string loading_dir = argv[1];

	//	Choose objective neuron in loop 1;
	int objective_neuron_index = atoi(argv[2]);
	string pre_neuron_filename = loading_dir + "preNeuron.txt";
	vector<int> pre_neuron_type;
	ReadColumn(pre_neuron_filename, 0, 6, pre_neuron_type);
	cout << ">> Target neuron in pre-network: ";
	if (pre_neuron_type[objective_neuron_index] == 1) {
		cout << "#" << objective_neuron_index << " neuron is an excitatory neuron." << endl;
	} else {
		cout << "#" << objective_neuron_index << " neuron is an inhibitory neuron." << endl;
	}

	string list_str =  argv[3]; 
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
	t_range[0] = atof(argv[4]);
	t_range[1] = atof(argv[5]);
	printf(">> Time range is (%.2f, %.2f] ms.\n", t_range[0], t_range[1]);

	string current_filename = loading_dir + "postI.txt";
	cout << ">> Calculating LFP ..." << endl;
	int total_neuron_number = atoi(argv[6]);
	vector<double> lfp;
	// LFP(t_range, total_neuron_number, list, potential_filename, excitatory_conductance_filename, inhibitory_conductance_filename, lfp);
	LFP(t_range, total_neuron_number, list, current_filename, lfp);

	//	Output data:
	string out_dir = "./lfp/file-txt/";
	string lfp_filename = out_dir + "lfp.txt";
	string raster_filename = out_dir + "raster.txt";

	//	Output LFP:
	cout << ">> Outputing LFP and spike train ..." << endl;
	OutputLFP(lfp, lfp_filename);
	//	Output raster of objective neruon:
	string original_raster_filename = loading_dir + "rasterPre.txt";
	vector<double> objective_neuron_raster;
	ReadLine(original_raster_filename, objective_neuron_index, objective_neuron_raster);
	OutputSpikeTrain(t_range, objective_neuron_raster, raster_filename);
	cout << ">> Finished." << endl;
	return 0;
}