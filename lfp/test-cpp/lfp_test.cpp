//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-08 21:33:10
//	Description: test file for lfp.h and lfp.cpp
//***************
#include "../lfp.h"
#include "../../io/io.h"
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

int main() {
	cout << "BEGIN: " << endl;
	//	Defined folder path;
	cout << ">> Enter existing target path for output: " << endl;
	string folder;
	cin >> folder;
	//	Choose objective neuron in loop 1;
	cout << ">> Preparing parameters..." << endl;
	int objective_neuron_index = 64;
	//	Choose objective time range;
	double t_range[2]; // t_range[0] = t_min; t_range[1] = t_max;
	t_range[0] = 1000;
	t_range[1] = 10000;
	//	Search in connectivity matrix and load connectivity vector of chosen neuron;
	cout << ">> Searching objective neuron..." << endl;
	string connectivity_path = folder + "conMat.txt";
	vector<int> neuronal_connectivity_list; // list of neuronal index from those which dirctly connected with objective neuron;
	Read1D(connectivity_path, objective_neuron_index, 0, neuronal_connectivity_list);
	//	Read neuron indices of neurons in the neuronal_connectivity_list;
	vector<int> connected_neurons; // indices of connected neurons;
	for (int i = 0; i < neuronal_connectivity_list.size(); i ++) {
		if (neuronal_connectivity_list[i] == 1) connected_neurons.push_back(i);
	}
	unique(connected_neurons.begin(),connected_neurons.end());
	sort(connected_neurons.begin(),connected_neurons.end());
	cout << connected_neurons.size() << endl;

	//	Search in files of neurons in loop 2, and load neuronal parameters, including potential, excitatory_conductance and inhibitory conductance, of neurons in the neuronal_connectivity_list;
	cout << ">> Loading parameters of connected neurons..." << endl;
	string potential_path = folder + "postV.txt";
	string excitatory_conductance_path = folder + "postGE.txt";
	string inhibitory_conductance_path = folder + "postGI.txt";

	vector<vector<double> > connected_neuron_potential, connected_neuron_excitatory_conductance, connected_neuron_inhibitory_conductance;
	cout << ">> Calculating LFP..." << endl;
	vector<double> lfp;
	LFP(t_range, connected_neurons, potential_path, excitatory_conductance_path, inhibitory_conductance_path, lfp);
	cout << lfp.size() << endl;
	//	Output data:
	string out_dir = "./";
	string lfp_path = out_dir + "lfp_test.txt";
	string raster_path = out_dir + "raster_test.txt";
	//	Output LFP:
	cout << ">> Outputing data..." << endl;
	OutLFP(lfp_path, lfp);
	//	Output raster of objective neruon:
	string original_raster_path = folder + "rasterPre.txt";
	vector<double> objective_neuron_raster;
	Read1D(original_raster_path, objective_neuron_index, 0, objective_neuron_raster);
	OuttSpikeTrain(raster_path, objective_neuron_raster, t_range);
	cout << "END!" << endl;
	return 0;
}
