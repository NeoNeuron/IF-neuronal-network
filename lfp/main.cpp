//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-05 14:44:33
//	Description: test file for lfp.h and lfp.cpp
//***************
#include "lfp.h"
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

int main() {
	cout << "BEGIN: " << endl;
	//	Defined folder path;
	cout << ">> Enter directory for neural data: ";
	string loading_dir;
	cin >> loading_dir;
	//	Choose objective neuron in loop 1;
	cout << ">> Preparing parameters..." << endl;
	cout << ">> >> Objective neuron index: ";
	int objective_neuron_index;
	cin >> objective_neuron_index;
	//	Choose objective time range;
	double t_range[2]; // t_range[0] = t_min; t_range[1] = t_max;
	cout << ">> >> Lower time limit (ms): ";
	cin >> t_range[0];
	cout << ">> >> Upper time limit (ms): ";
	cin >> t_range[1];
	int connecting_order;
	cout << ">> >> Connecting order: ";
	cin >> connecting_order;
	//	Search in connectivity matrix and load connectivity vector of chosen neuron;
	cout << ">> Searching objective neuron ..." << endl;
	string connectivity_filename_name = loading_dir + "conMat.txt";
	vector<int> neuronal_connectivity_list; // list of neuronal index from those which dirctly connected with objective neuron;
	ReadLine(connectivity_filename_name, objective_neuron_index, neuronal_connectivity_list);
	//	Read neuron indices of neurons in the neuronal_connectivity_list;
	vector<int> connected_neurons;
	for (int i = 0; i < neuronal_connectivity_list.size(); i ++) {
		if (neuronal_connectivity_list[i] == 1) connected_neurons.push_back(i);
	}
	int first_order_size = connected_neurons.size();
	if (connecting_order == 2) {
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
	}
	cout << ">> The number of connected neurons is " << connected_neurons.size() << "." << endl;

	int classification;
	cout << ">> Classification of sub neuron cluster: " << endl << "[1]Excitatory; [2]Inhibitory; [3]Both; \t";
	cin >> classification;
	vector<int> post_neuron_type;
	string neuron_filename = loading_dir + "postNeuron.txt";
	ReadColumn(neuron_filename, 0, 6, post_neuron_type);
	int initial_size = connected_neurons.size();
	vector<int> connected_neurons_copy = connected_neurons;
	switch(classification) {
		case 1: {
			connected_neurons.clear();
			for (int i = 0; i < initial_size; i++) { // loop for seleted neurons in first step;
				if (post_neuron_type[connected_neurons_copy[i]] == 1) connected_neurons.push_back(connected_neurons_copy[i]);
			}
			break;
		}
		case 2: {
			connected_neurons.clear();
			for (int i = 0; i < initial_size; i++) { // loop for seleted neurons in first step;
				if (post_neuron_type[connected_neurons_copy[i]] == 0) connected_neurons.push_back(connected_neurons_copy[i]);
			}
		}
		case 3: break;
		default: break;
	}
	connected_neurons_copy.clear();
	
	//	Search in files of neurons in loop 2, and load neuronal parameters, including potential, excitatory_conductance and inhibitory conductance, of neurons in the neuronal_connectivity_list;
	cout << ">> Loading parameters of connected neurons ..." << endl;
	string potential_filename = loading_dir + "postV.txt";
	string excitatory_conductance_filename = loading_dir + "postGE.txt";
	string inhibitory_conductance_filename = loading_dir + "postGI.txt";
	
	cout << ">> Calculating LFP ..." << endl;
	vector<double> lfp;
	LFP(t_range, connected_neurons, potential_filename, excitatory_conductance_filename, inhibitory_conductance_filename, lfp);

	//	Output data:
	string out_dir = "./file_txt/";
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
	cout << ">> END!" << endl;
	return 0;
}