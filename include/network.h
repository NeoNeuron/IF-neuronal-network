//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2018-05-31
//	Description: define Struct SpikeElement and Class NeuronalNetwork;
//***************
#ifndef _IFNET_NETWORK_H_
#define _IFNET_NETWORK_H_

#include "io.h"
#include "neuron.h"
#include "get-config.h"
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
	int neuron_number_;	// number of the neurons in the group;
	vector<double*> dym_vals_; // dynamic variables of neurons;
	vector<double*> dym_vals_new_; // temporal dynamic variables of neurons;
	vector<Neuron> neurons_;
	vector<bool> types_; // vector to store types of neurons;
	vector<vector<bool> > con_mat_; // built-in matrix for neuronal connectivity;
	bool is_con_;
	vector<vector<double> > external_exc_inputs_; // temp storage of external Poisson input;
	vector<vector<double> > external_inh_inputs_;
	double interaction_delay_;
	// Functions:
	// Sort spikes within single time interval, and return the time of first spike;
	double SortSpikes(vector<double*> &dym_vals_new, vector<int>& update_list, vector<int>& fired_list, double t, double dt, vector<SpikeElement> &T);

public:
	//	Neuronal network initialization:
	NeuronalNetwork(int neuron_number) {
		neuron_number_ = neuron_number;
		dym_vals_.resize(neuron_number_);
		dym_vals_new_.resize(neuron_number_);
		for (int i = 0; i < neuron_number_; i++) {
			neurons_.push_back(Neuron(dym_vals_[i]));
			//TODO: allocate memory for dynamic variable for different neuronal model;
			//dym_vals_new[i] = new double[GetLen(dym_vals_[i])];
			dym_vals_new_[i] = new double[4];
		}
		types_.resize(neuron_number_, false);
		con_mat_.resize(neuron_number_, vector<bool>(neuron_number_, false));
		is_con_ = false;
		interaction_delay_ = 0.0;
	}
	
	~NeuronalNetwork() {
		for (int i = 0; i < neuron_number_; i++) {
			delete [] dym_vals_new_[i];
		}
	}
	// Initialize network connectivity matrix:
	// Three options:
	// 0. given pre-defined network connectivity matrix;
	// 1. small-world network, defined by connectivity density and rewiring;
	// 2. randomly connected network;
	void InitializeConnectivity(map<string, string> &m_config, string prefix);
	
	// INPUTS:
	// Set interneuronal coupling strength;
	void SetS(bool function, double val);

	// Set feedforward inputing strength: (default: 5e-3)
	void SetF(bool function, double val);

	// Set time period of refractory:
	void SetRef(double t_ref);

  // Set interaction delay between neurons;
	void SetDelay(double val) { interaction_delay_ = val; }

	// 	Initialize neuronal types in the network;
	//	p: the probability of the presence of excitatory neuron;
	//	seed: seed for random number generator;
	void InitializeNeuronalType(double p, int seed);

	//	Set driving type: true for external Poisson driven, false for internal ones;
	void SetDrivingType(bool driving_type);

	//	Initialize internal homogeneous feedforward Poisson rate;
	//	rates: Poisson rates; Each line represents the two parameters for each neuron,
	//		the first if excitatory	diring rate, and the second is inhibitory.
	void InitializeInternalPoissonRate(vector<vector<double> >& rates);

	//	Initialize external homogeneous feedforward Poisson process;
	//	rates: Poisson rates; Each line represents the two parameters for each neuron,
	//		the first if excitatory	diring rate, and the second is inhibitory.
	//	tmax: maximum time range for Poisson process;
	//	seed: seed for built-in random generator;
	void InitializeExternalPoissonProcess(vector<vector<double> >& rates, double tmax, int seed);

	// 	Input new spikes for neurons all together;
	void InNewSpikes(vector<vector<Spike> > &data);

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

	// Print cycle:
	void PrintCycle();

	//	Output potential to *.csv file;
	void OutPotential(FILEWRITE& file);

  //	Output synaptic conductance to *.csv files:
  //		BOOL function: function of synaptic conductance, true for excitation, false for inhibition;
  void OutConductance(FILEWRITE& file, bool function);

  //	Output the total membrane ionic current of each neuron:
  void OutCurrent(FILEWRITE& file);

  //	Output the partial membrane ionic current of each neuron:
  void OutPartialCurrent(FILEWRITE& file, bool type);
	
	// Save connectivity matrix
	void SaveConMat(string connecting_matrix_file);

	//	Save corrent neuronal States:
	//	Define a 'neuronFile' type to store neuronal condition;
	//	A ROW VECTOR:
	//		0: neuronal type;
	//		1: neuronal index;
	//		2: membrane potential;
	//		3: excitatory conductivity;
	//		4: inhibitory conductivity;
	//		5: remaining refractory period;
	void SaveNeuron(string neuron_file);

	int OutSpikeTrains(string path);

  //  Output spikes before t, including their modes and functions;
	void GetNewSpikes(double t, vector<vector<Spike> >& data);

	void GetNeuronType(vector<bool> & x);

	int GetNeuronNumber();

	void GetConductance(int i, bool function);
};

void UpdateSystemState(NeuronalNetwork & pre_network, NeuronalNetwork & post_network, vector<vector<bool> > &connectivity_matrix, double t, double dt);

#endif // _IFNET_NETWORK_H_
