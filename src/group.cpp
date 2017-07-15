//******************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Description: Define class Neuron, structure Spike and NeuronState;
//	Date: 2017-02-21 16:06:30
//******************************
#include "../include/group.h"
#include "../include/io.h"
#include <iostream>
#include <algorithm>
#include <ctime>
#include <cmath>

using namespace std;

bool Compare(const SpikeElement &x, const SpikeElement &y) {
	return x.t < y.t;
}

double NeuronalNetwork::SortSpikes(double t, double dt, vector<SpikeElement>& T) {
	double SET = -1;
	SpikeElement ADD;
	for (int i = 0; i < neuron_number_; i++) {
		/*if (i == 0) {
			_SET = neurons_[i].TemporallyUpdateNeuronalState(t, dt, true);
		}*/
		SET = neurons_[i].TemporallyUpdateNeuronalState(t, dt, external_excitatory_inputs_[i], external_inhibitory_inputs_[i]);
		if (SET > 0) {
			ADD.index = neurons_[i].GetNeuronIndex();
			ADD.t = SET;
			ADD.type = neurons_[i].GetNeuronalType();
			T.push_back(ADD);
		}
	}
	if (T.size() == 0) {
		return -1;
	} else if (T.size() == 1) {
		return (T.front()).t;
	} else {
		sort(T.begin(), T.end(), Compare);
		return (T.front()).t;
	}
}

void NeuronalNetwork::SetConnectingDensity(int density){
	connecting_density_ = density;
	connectivity_matrix_.SetConnectingDensity(connecting_density_);
}

void NeuronalNetwork::InitializeNeuronalType(double p, int seed) {
	srand(seed);
	double x = 0;
	int counter = 0;
	for (int i = 0; i < neuron_number_; i++) {
		x = rand() / (RAND_MAX + 1.0);
		if (x < p) {
			neurons_[i].SetNeuronType(true);
			counter++;
		} else neurons_[i].SetNeuronType(false);
	}
	printf(">> There are %d excitatory neurons and %d inhibitory neurons ", counter, neuron_number_-counter);
}

void NeuronalNetwork::InitializeInternalPoissonRate(bool function, double rate) {
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].SetPoissonRate(function, rate);
	}
}

void NeuronalNetwork::InitializeExternalPoissonProcess(bool function, double rate_excitatory, double rate_inhibitory, double tmax, int seed) {
	vector<double> ADD;
	ADD.push_back(0);
	srand(seed);
	double x; // temp random number;
	double tLast; // last existing data point;
	if (function == true) {
		//ofstream external_poisson;
		//external_poisson.open("external_poisson_new.txt");
		external_excitatory_inputs_.resize(neuron_number_, ADD);
		// for (int i = 0; i < neuron_number_; i++) {
		// 	external_excitatory_inputs_.push_back(ADD);
		// }
		double rate;
		for (int i = 0; i < neuron_number_; i++) {
			tLast = external_excitatory_inputs_[i].front();
			if (neurons_[i].GetNeuronalType() == true) rate = rate_excitatory;
			else rate = rate_inhibitory;
			while (tLast < tmax) {
				x = rand() / (RAND_MAX + 1.0);
				while (log(x) < -1e12) x = rand() / (RAND_MAX + 1.0);
				tLast -= log(x) / rate;
				external_excitatory_inputs_[i].push_back(tLast);
				//external_poisson << (double)tLast << ",";
			}
			//external_poisson << endl;
		}
	} else {
		external_inhibitory_inputs_.resize(neuron_number_, ADD);
		// for (int i = 0; i < neuron_number_; i++) {
		// 	external_inhibitory_inputs_.push_back(ADD);
		// }
		double rate;
		for (int i = 0; i < neuron_number_; i++) {
			if (neurons_[i].GetNeuronalType() == true) rate = rate_excitatory;
			else rate = rate_inhibitory;
			tLast = external_inhibitory_inputs_[i].front();
			while (tLast < tmax) {
				x = rand() / (RAND_MAX + 1.0);
				tLast -= log(x) / rate;
				external_inhibitory_inputs_[i].push_back(tLast);
			}
		}
	}
}

void NeuronalNetwork::InNewSpikes(vector<vector<Spike> > & data) {
	for (int i = 0; i < neuron_number_; i++) {
		if (data[i].size() != 0) {
			for (vector<Spike>::iterator it = data[i].begin(); it != data[i].end(); it++) {
				neurons_[i].InSpike(*it);
			}
		}
	}
}

void NeuronalNetwork::LoadNeuronalState(string neuron_file) {
	vector<vector<double> > neuronalSetups;
	Read2D(neuron_file, neuronalSetups);
	NeuronalState add;
	for (int i = 0; i < neuron_number_; i++) {
		if (neuronalSetups[i][0] == 1) add.type = true;
		else add.type = false;
		add.index = neuronalSetups[i][1];
		add.membrane_potential_ = neuronalSetups[i][2];
		add.ge = neuronalSetups[i][3];
		add.gi = neuronalSetups[i][4];
		add.remaining_refractory_time = neuronalSetups[i][5];
		neurons_[i].LoadNeuronalState(add);
	}
}

void NeuronalNetwork::LoadConnectivityMat(string connecting_matrix_file) {
	vector<vector<int> > connecting_matrix;
	Read2D(connecting_matrix_file, connecting_matrix);
	connectivity_matrix_.LoadMatrix(connecting_matrix);
}

void NeuronalNetwork::SetDrivingType(bool driving_type) {
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].SetDrivingType(driving_type);
	}
}

void NeuronalNetwork::Rewire(double p, int seed, bool output_option) {
	connectivity_matrix_.Rewire(p, seed, output_option);
	if (output_option == true) {
		double mean_path, mean_clustering_coefficient;
		mean_path = connectivity_matrix_.GetMeanPath();
		mean_clustering_coefficient = connectivity_matrix_.GetMeanClusteringCoefficient();
		cout << "After rewiring," << endl;
		cout << "The mean characteristic path is " << setprecision(4) << (double)mean_path << "." << endl;
		cout << "The mean clustering coefficient is " << setprecision(4) << (double)mean_clustering_coefficient << "." << endl;
	}
}

void NeuronalNetwork::UpdateNetworkState(double t, double dt) {
	vector<SpikeElement> T;
	double newt;
	newt = SortSpikes(t, dt, T);
	if (newt < 0) {
		for (int i = 0; i < neuron_number_; i++) {
			neurons_[i].UpdateNeuronalState(t, dt, external_excitatory_inputs_[i], external_inhibitory_inputs_[i]);
		}
	} else {
		while (newt > 0) {
			int IND = (T.front()).index;
			Spike ADD_mutual;
			ADD_mutual.mode = false;
			ADD_mutual.function = (T.front()).type;
			ADD_mutual.t = newt;
			for (int j = 0; j < neuron_number_; j++) {
				if (j == IND) {
					neurons_[j].Fire(t, newt - t);
				} else {
					neurons_[j].UpdateNeuronalState(t, newt - t, external_excitatory_inputs_[j], external_inhibitory_inputs_[j]);
					if (connectivity_matrix_.ReadMatrix(IND,j) != 0) {
						neurons_[j].InSpike(ADD_mutual);
					}
				}
			}
			dt = t + dt - newt;
			t = newt;
			T.clear();
			newt = SortSpikes(t, dt, T);
		}
		for (int i = 0; i < neuron_number_; i++) {
			neurons_[i].UpdateNeuronalState(t, dt, external_excitatory_inputs_[i], external_inhibitory_inputs_[i]);
		}
	}
}

void NeuronalNetwork::OutPotential(string path) {
	vector<double> potential(neuron_number_);
	for (int i = 0; i < neuron_number_; i++) {
		potential[i] = neurons_[i].GetPotential();
	}
	Print1D(path, "app", 0, potential);
}

void NeuronalNetwork::OutConductance(string path, bool function) {
	vector<double> conductance(neuron_number_);
	if (function) {
		for (int i = 0; i < neuron_number_; i++) {
			conductance[i] = neurons_[i].GetConductance(true);
		}
	} else {
		for (int i = 0; i < neuron_number_; i++) {
			conductance[i] = neurons_[i].GetConductance(false);
		}
	}
  Print1D(path, "app", 0, conductance);
}

void NeuronalNetwork::OutCurrent(string path) {
	vector<double> current(neuron_number_);
	for (int i = 0; i < neuron_number_; i++) {
		current[i] = neurons_[i].OutTotalCurrent();
	}
	Print1D(path, "app", 0, current);
}

void NeuronalNetwork::OutPartialCurrent(string path, bool type) {
	vector<double> current(neuron_number_);
	for (int i = 0; i < neuron_number_; i++) {
		current[i] = neurons_[i].OutLeakyCurrent() + neurons_[i].OutSynapticCurrent(type);
	}
	Print1D(path, "app", 0, current);
}


void NeuronalNetwork::Save(string neuron_file, string connecting_matrix_file) {
	ofstream data;
	data.open(neuron_file.c_str());
	NeuronalState add;
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].Save(add);
		data << (bool)add.type << ',';
		data << (int)add.index << ',';
		data << (double)add.membrane_potential_ << ',';
		data << (double)add.ge << ',';
		data << (double)add.gi << ',';
		data << (double)add.remaining_refractory_time << endl;
	}
	data.close();
	connectivity_matrix_.OutMatrix(connecting_matrix_file);
}


void NeuronalNetwork::OutSpikeTrains(string path) {
	vector<vector<double> > spikes(neuron_number_);
	vector<double> add_spike_train;
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].OutSpikeTrain(add_spike_train);
		spikes[i] = add_spike_train;
	}
	Print2D(path, "trunc", spikes);
}

void NeuronalNetwork::GetNewSpikes(double t, vector<vector<Spike> >& data) {
	data.clear();
	data.resize(neuron_number_);
	vector<Spike> x;
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].GetNewSpikes(t, x);
		data[i] = x;
	}
}

void NeuronalNetwork::GetNeuronType(vector<bool>& types) {
	types.clear();
	types.resize(neuron_number_);
	for (int i = 0; i < neuron_number_; i++) {
		types[i] = neurons_[i].GetNeuronalType();
	}
}

int NeuronalNetwork::GetNeuronNumber() {
	return neuron_number_;
}

void NeuronalNetwork::GetConductance(int i, bool function) {
	neurons_[i].GetConductance(function);
}

void NeuronalNetwork::RestoreNeurons() {
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].Reset();
	}
	external_excitatory_inputs_.clear();
	external_inhibitory_inputs_.clear();
}

void UpdateSystemState(NeuronalNetwork & pre_network, NeuronalNetwork & post_network, vector<vector<bool> > & connectivity_matrix, double t, double dt) {
	//	Update pre_network in two network system;
	pre_network.UpdateNetworkState(t, dt);
	//	Transmit spiking information from pre_network to post_network;
	vector<vector<Spike> > tempPreSpikes, tempPostSpikes;
	pre_network.GetNewSpikes(t, tempPreSpikes);
	// Set transmit time:
	// double transmit_time = 2;
	// for (vector<vector<Spike> >::iterator it = tempPreSpikes.begin(); it != tempPreSpikes.end(); it ++) {
	// 	if (it -> size() > 0) {
	// 		for (vector<Spike>::iterator itt = it -> begin(); itt != it -> end(); itt ++) {
	// 			itt -> t += transmit_time;
	// 		}
	// 	}
	// }
	// find temporal spiking sequence for post network;
	tempPostSpikes.resize(connectivity_matrix[1].size());
	for (int i = 0; i < connectivity_matrix.size(); i++) {
		if (tempPreSpikes[i].size() != 0) {
			for (int j = 0; j < connectivity_matrix.size(); j++) {
				if (connectivity_matrix[i][j] == true) {
					tempPostSpikes[j].insert(tempPostSpikes[j].end(), tempPreSpikes[i].begin(), tempPreSpikes[i].end());
				}
			}
		}
	}
	post_network.InNewSpikes(tempPostSpikes);
	//	Update post_network in two network system;
	post_network.UpdateNetworkState(t, dt);
}
