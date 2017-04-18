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

using namespace std;

//	Function of calculating LFP with point current source model in 1-D loop network case;
//	arguments:
//	argv[1] = loading directory for neural data;
//	argv[2] = neuronal index for target neuron in pre network;
//	argv[3] = connecting order of neurons generating LFP;
//	argv[4] = lower bond of time range;
//	argv[5] = upper bond of time range;
//	argv[6] = detailed connecting classification of neurons in post network;
//		"all" = all neuron in given order;
//		"excitatory" = excitatory neurons among neurons in given order;
//		"inhibitory" = inhibitory neurons among neurons in given order;
//	argv[7] = number of neurons that chosen in the subset of given classification; if argv[7] = 0, all neurons that pass previous selections are preserved;
//	argv[8] = total neuron number;
int main(int argc, const char* argv[]) {
	if (argc != 9) {
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
	int connecting_order =  atoi(argv[3]);
	printf(">> Connecting order is %d.\n", connecting_order);

	//	Choose objective time range;
	double t_range[2]; // t_range[0] = t_min; t_range[1] = t_max;
	t_range[0] = atof(argv[4]);
	t_range[1] = atof(argv[5]);
	printf(">> Time range is (%.2f, %.2f] ms.\n", t_range[0], t_range[1]);
	//	Search in connectivity matrix and load connectivity vector of chosen neuron;
	string connectivity_filename_name = loading_dir + "conMat.txt";
	vector<int> neuronal_connectivity_list; // list of neuronal index from those which dirctly connected with objective neuron;
	ReadLine(connectivity_filename_name, objective_neuron_index, neuronal_connectivity_list);
	//	Read neuron indices of neurons in the neuronal_connectivity_list;
	//	First order neuron indices;
	vector<int> connected_neurons;
	for (int i = 0; i < neuronal_connectivity_list.size(); i ++) {
		if (neuronal_connectivity_list[i] == 1) connected_neurons.push_back(i);
	}
	int first_order_size = connected_neurons.size();
	//	Second order neuron indices;
	if (connecting_order == 2) {
		vector<int> first_order_connected_neurons = connected_neurons;
		vector<int> second_order_connected_neurons;
		for (int i = 0; i < first_order_size; i ++) {
			neuronal_connectivity_list.clear();
			ReadLine(connectivity_filename_name, connected_neurons[i], neuronal_connectivity_list);
			for (int j = 0; j < neuronal_connectivity_list.size(); j ++) {
				if (neuronal_connectivity_list[j] == 1) second_order_connected_neurons.push_back(j);
			}
		}
		connected_neurons.clear();
		sort(second_order_connected_neurons.begin(),second_order_connected_neurons.end());
		unique_copy(second_order_connected_neurons.begin(),second_order_connected_neurons.end(), back_inserter(connected_neurons));
		second_order_connected_neurons.clear();
		second_order_connected_neurons = connected_neurons;
		connected_neurons.clear();
		set_difference(second_order_connected_neurons.begin(), second_order_connected_neurons.end(), first_order_connected_neurons.begin(), first_order_connected_neurons.end(), inserter(connected_neurons, connected_neurons.begin()));
	}
	int neuron_num = connected_neurons.size();
	printf(">> %d neurons in post-network connect(s) with target neuron in given order.\n", neuron_num);

	string classification = argv[6];
	cout << ">> Classification of sub neuron cluster is " << classification << endl;
	int num = atoi(argv[7]);
	vector<int> post_neuron_type;
	string post_neuron_filename = loading_dir + "postNeuron.txt";
	ReadColumn(post_neuron_filename, 0, 6, post_neuron_type);
	vector<neuron_type> types(neuron_num);
	for (int i = 0; i < neuron_num; i++) {
		types[i].index = connected_neurons[i];
		if (post_neuron_type[types[i].index] == 1) {
			types[i].type = true;
		} else {
			types[i].type = false;
		}
	}
	KeySelect(classification, types, connected_neurons);
	vector<int> temp_neurons;
	Sample(connected_neurons, temp_neurons, num);
	connected_neurons.clear();
	connected_neurons = temp_neurons;	
	// int target = connected_neurons[num];
	// connected_neurons.clear();
	// connected_neurons.push_back(target);
	//connected_neurons.erase(connected_neurons.begin() + num, connected_neurons.begin() + num + 1);
	// for (int i = 0; i < num; i ++) cout << connected_neurons[i] << '\t';
	// cout << endl; 
	cout << ">> The number of chosen neurons is " << connected_neurons.size() << "." << endl;
	sort(connected_neurons.begin(), connected_neurons.end());
	neuron_num = connected_neurons.size();
	printf(">> %d neurons in post-network is(are) chosen to generate local field potential.\n", neuron_num);
	
	//	Search in files of neurons in loop 2, and load neuronal parameters, including potential, excitatory_conductance and inhibitory conductance, of neurons in the neuronal_connectivity_list;
	cout << ">> Loading parameters of connected neurons ..." << endl;
	// string potential_filename = loading_dir + "postV.txt";
	// string excitatory_conductance_filename = loading_dir + "postGE.txt";
	// string inhibitory_conductance_filename = loading_dir + "postGI.txt";
	string current_filename = loading_dir + "postI.txt";
	  
	cout << ">> Calculating LFP ..." << endl;
	int total_neuron_number = atoi(argv[8]);
	vector<double> lfp;
	// LFP(t_range, total_neuron_number, connected_neurons, potential_filename, excitatory_conductance_filename, inhibitory_conductance_filename, lfp);
	LFP(t_range, total_neuron_number, connected_neurons, current_filename, lfp);

	//	Output data:
	string out_dir = "./lfp/file-txt/";
	string lfp_filename = out_dir + "lfp_test.txt";
	string raster_filename = out_dir + "raster_test.txt";

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