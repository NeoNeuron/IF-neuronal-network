//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2018-05-31
//	Description: define Struct SpikeElement and Class NeuronalNetwork;
//***************
#ifndef _IFNET_NETWORK_H_
#define _IFNET_NETWORK_H_

#include "neuron.h"
#include "get-config.h"
#include "connectivity_matrix.h"
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

struct SpikeElement {
	int index;	// The sequence order of spikes within single time interval;
	double t;	// exact spiking time;
	bool type;	// The type of neuron that fired;
};

class NeuronalNetwork {
private:
	// Parameters:
	Neuron *neurons_;
	int neuron_number_;	// number of the neurons in the group;
	ConnectivityMatrix connectivity_matrix_;
	vector<vector<double> > external_excitatory_inputs_; // temp storage of external Poisson input;
	vector<vector<double> > external_inhibitory_inputs_;
	double interaction_delay_;
	int connecting_density_;
	// Functions:
	// Sort spikes within single time interval, and return the time of first spike;
	double SortSpikes(double t, double dt, vector<SpikeElement> &T);

public:
	//	Neuronal network initialization:
	NeuronalNetwork(int neuron_number) {
		neuron_number_ = neuron_number;
		neurons_ = new Neuron[neuron_number_];
		for (int i = 0; i < neuron_number_; i++) neurons_[i].SetNeuronIndex(i);
		connectivity_matrix_.SetNeuronNumber(neuron_number_);
		interaction_delay_ = 0.0;
		connecting_density_ = 0;
	}
	
	// Initialize network connectivity matrix:
	// Three options:
	// 0. given pre-defined network connectivity matrix;
	// 1. small-world network, defined by connectivity density and rewiring;
	// 2. randomly connected network;
	void InitializeConnectivity(map<string, string> &m_config);
	
	void RandNet(double p, int seed) {
		connectivity_matrix_.RandNet(p, seed);
	}
	// INPUTS:
	// Set interneuronal coupling strength;
	void SetS(bool function, double val);

  // Set interaction delay between neurons;
	void SetDelay(double val);

	// 	Initialize neuronal types in the network;
	//	DOUBLE p: the probability of the presence of excitatory neuron;
	//	DOUBLE seed: seed for random number generator;
	void InitializeNeuronalType(double p, int seed);

	//	Set driving type: true for external Poisson driven, false for internal ones;
	void SetDrivingType(bool driving_type);

	//	Initialize internal homogeneous feedforward Poisson rate;
	//	DOUBLE rates: Poisson rates; Each line represents the two parameters for each neuron,
	//		the first if excitatory	diring rate, and the second is inhibitory.
	void InitializeInternalPoissonRate(vector<vector<double> >& rates);

	//	Initialize external homogeneous feedforward Poisson process;
	//	DOUBLE rates: Poisson rates; Each line represents the two parameters for each neuron,
	//		the first if excitatory	diring rate, and the second is inhibitory.
	//	DOUBLE tmax: maximum time range for Poisson process;
	//	INT seed: seed for built-in random generator;
	void InitializeExternalPoissonProcess(vector<vector<double> >& rates, double tmax, int seed);

	// Set feedforward inputing strength: (default: 5e-3)
	void SetF(bool function, double val);

	// 	Input new spikes for neurons all together;
	void InNewSpikes(vector<vector<Spike> > &data);

	//Load previous setups;

	// Load neuronal states:
	// Define a 'NeuronState' type to store neuronal condition;
	// A ROW VECTOR:
	//	0: neuronal type;
	//	1: neuronal index;
	//	2: membrane potential;
	//	3: excitatory conductivity;
	//	4: inhibitory conductivity;
	//	5: remaining refractory period;
	void LoadNeuronalState(string neuronFile);

	// DYNAMICS:

	//	Restore neuronal state for all neurons, including neuronal potential, conductances, refractory periods and external network drive;
	void RestoreNeurons();

	//	Update network state:
	void UpdateNetworkState(double t, double dt);

	// OUTPUTS:

	//	Output potential to *.csv file;
	void OutPotential(string path);

  //	Output synaptic conductance to *.csv files:
  //		BOOL function: function of synaptic conductance, true for excitation, false for inhibition;
  void OutConductance(string path, bool function);

  //	Output the total membrane ionic current of each neuron:
  void OutCurrent(string path);

  //	Output the partial membrane ionic current of each neuron:
  void OutPartialCurrent(string path, bool type);

	//	Save corrent neuronal States and connectivity matrix:
	//	Define a 'neuronFile' type to store neuronal condition;
	//	A ROW VECTOR:
	//		0: neuronal type;
	//		1: neuronal index;
	//		2: membrane potential;
	//		3: excitatory conductivity;
	//		4: inhibitory conductivity;
	//		5: remaining refractory period;
	void Save(string neuron_file, string connecting_matrix_file);

	int OutSpikeTrains(string path);

  //  Output spikes before t, including their modes and functions;
	void GetNewSpikes(double t, vector<vector<Spike> >& data);

	void GetNeuronType(vector<bool> & x);

	int GetNeuronNumber();

	void GetConductance(int i, bool function);
};

void UpdateSystemState(NeuronalNetwork & pre_network, NeuronalNetwork & post_network, vector<vector<bool> > &connectivity_matrix, double t, double dt);

#endif // _IFNET_NETWORK_H_
