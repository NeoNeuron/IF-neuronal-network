//******************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Description: Define class Neuron, structure Spike and NeuronState;
//	Date: 2017-02-21 16:06:30
//******************************
#include "../include/network.h"
#include <iostream>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <cstring>

using namespace std;

bool Compare(const SpikeElement &x, const SpikeElement &y) {
	return x.t < y.t;
}

void Scan(vector<bool> & mat, int target_value, vector<int> &output_indices) {
	output_indices.clear();
	for (int s = 0; s < mat.size(); s++) {
		if (mat[s] == target_value) output_indices.push_back(s);
	}
}

void NeuronalNetwork::InitializeConnectivity(map<string, string> &m_config, string prefix) {
	int connecting_mode = atoi(m_config[prefix + "ConnectingMode"].c_str());
	if (connecting_mode == 0) { // External connectivity matrix;
		vector<vector<int> > connecting_matrix;
		Read2D(m_config[prefix + "MatPath"], connecting_matrix);
		if (connecting_matrix.size() != neuron_number_ || connecting_matrix[0].size() != neuron_number_) {
			throw runtime_error("wrong size of connectivity matrix");
		} else {
			for (size_t i = 0; i < neuron_number_; i ++) {
				for (size_t j = 0; j < neuron_number_; j ++) {
					if (connecting_matrix[i][j]) {
						con_mat_[i][j] = true;
						if (!is_con_) is_con_ = true;
					}
				}
			}
		}
	} else if (connecting_mode == 1) {
		int con_density = atoi(m_config[prefix + "ConnectingDensity"].c_str());
		for (int i = 0; i < neuron_number_; i++)  {
			for (int j = 0; j < neuron_number_; j++) {
				if (i != j) {
					if (abs(i - j) <= con_density or neuron_number_ - abs(i - j) <= con_density) {
					con_mat_[i][j] = true;
					if (!is_con_) is_con_ = true;
					}
				}
			}
		}
		double rewiring_probability = atof(m_config[prefix + "RewiringProbability"].c_str());
		bool output_option;
		istringstream(m_config[prefix + "PrintRewireResult"]) >> boolalpha >> output_option;
		// Generate networks;
		cout << 2 * neuron_number_ * con_density << " connections total with ";
		srand(atoi(m_config[prefix + "NetSeed"].c_str()));
		double x; // random variable;
		int ind, empty_connection, count = 0;
		vector<int> ones, zeros;
		for (int i = 0; i < neuron_number_; i++) {
			Scan(con_mat_[i], 1, ones);
			for (int j = 0; j < ones.size(); j++) {
				x = rand() / (RAND_MAX + 1.0);
				if (x <= rewiring_probability) {
					Scan(con_mat_[i], 0, zeros);
					for (vector<int>::iterator it = zeros.begin(); it != zeros.end(); it++) {
						if (*it == i) {
							zeros.erase(it);
							break;
						}
					}
					empty_connection = zeros.size();
					ind = rand() % empty_connection;
					con_mat_[i][zeros[ind]] = true;
					con_mat_[i][ones[j]] = false;
					count += 1;
				}
			}
		}
		cout << count << " rewirings." << endl;
	} else if (connecting_mode == 2) {
		double con_prob_ee = atof(m_config[prefix + "ConnectingProbabilityEE"].c_str());
		double con_prob_ei = atof(m_config[prefix + "ConnectingProbabilityEI"].c_str());
		double con_prob_ie = atof(m_config[prefix + "ConnectingProbabilityIE"].c_str());
		double con_prob_ii = atof(m_config[prefix + "ConnectingProbabilityII"].c_str());
		srand(atoi(m_config[prefix + "NetSeed"].c_str()));
		double x;
		for (size_t i = 0; i < neuron_number_; i ++) {
			for (size_t j = 0; j < neuron_number_; j ++) {
				x = rand() / (RAND_MAX * 1.0);
				if (types_[i]) {
					if (types_[j]) {
						if (x <= con_prob_ee) con_mat_[i][j] = true;
					} else {
						if (x <= con_prob_ei) con_mat_[i][j] = true;
					}
				} else {
					if (types_[j]) {
						if (x <= con_prob_ie) con_mat_[i][j] = true;
					} else {
						if (x <= con_prob_ii) con_mat_[i][j] = true;
					}
				}
				if (con_mat_[i][j]) {
					if (!is_con_) is_con_ = true;
				}
			}
		}
	}
}

bool CheckExist(int index, vector<int> &list) {
	for (int i = 0; i < list.size(); i ++) {
		if (list[i] == index) return true;
	}
	return false;
}

double NeuronalNetwork::SortSpikes(vector<double*> &dym_vals_new, vector<int> &update_list, vector<int> &fired_list, double t, double dt, vector<SpikeElement> &T) {
	double SET;
	SpikeElement ADD;	
	// start scanning;
	double id;
	for (int i = 0; i < update_list.size(); i++) {
		id = update_list[i];
		// Check whether id's neuron is in the fired list;
		if (CheckExist(id, fired_list)) {
			for (int j = 1; j < 3; j++) dym_vals_new[id][j] = dym_vals_[id][j];
			neurons_[id].UpdateConductance(dym_vals_new[id], t, dt);
		} else {
			SET = neurons_[id].TemporallyUpdateNeuronalState(dym_vals_[id], dym_vals_new[id], t, dt, external_exc_inputs_[id], external_inh_inputs_[id]);
			if (SET >= 0) {
				ADD.index = id;
				ADD.t = SET;
				ADD.type = neurons_[id].GetNeuronalType();
				T.push_back(ADD);
			}
		}
	}
	if (T.empty()) {
		return -1;
	} else if (T.size() == 1) {
		return (T.front()).t;
	} else {
		sort(T.begin(), T.end(), Compare);
		return (T.front()).t;
	}
}

void NeuronalNetwork::SetS(bool function, double val) {
	for (int i = 0; i < neuron_number_; i ++) {
		neurons_[i].SetSynapticStrength(function, val);
	}
}

void NeuronalNetwork::SetF(bool function, double val) {
	for (int i = 0; i < neuron_number_; i ++) {
		neurons_[i].SetFeedforwardStrength(function, val);
	}
}

void NeuronalNetwork::SetRef(double t_ref) {
	for (int i = 0; i < neuron_number_; i ++) { neurons_[i].SetRef(t_ref); }
}

void NeuronalNetwork::InitializeNeuronalType(double p, int seed) {
	srand(seed);
	double x = 0;
	int counter = 0;
	for (int i = 0; i < neuron_number_; i++) {
		x = rand() / (RAND_MAX + 1.0);
		if (x < p) {
			neurons_[i].SetNeuronType(true);
			types_[i] = true;
			counter++;
		} else neurons_[i].SetNeuronType(false);
	}
	printf(">> %d excitatory and %d inhibitory neurons ", counter, neuron_number_-counter);
}

void NeuronalNetwork::InitializeInternalPoissonRate(vector<vector<double> >& rates) {
	if (rates.size() != neuron_number_) {
		cout << "Error inputing length! (Not equal to the number of neurons in the net)";
		return;
	}
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].SetPoissonRate(true, rates[i][0]);
		neurons_[i].SetPoissonRate(false, rates[i][1]);
	}
}

void NeuronalNetwork::InitializeExternalPoissonProcess(vector<vector<double> >& rates, double tmax, int seed) {
	if (rates.size() != neuron_number_) {
		cout << "Error inputing length! (Not equal to the number of neurons in the net)";
		return;
	}
	vector<double> ADD = {0.0};
	srand(seed);
	double x; // temp random number;
	double tLast; // last existing data point;
	double rate;
	external_exc_inputs_.resize(neuron_number_, ADD);
	external_inh_inputs_.resize(neuron_number_, ADD);
	for (int i = 0; i < neuron_number_; i++) {
		rate = rates[i][0];
		tLast = external_exc_inputs_[i].front();
		while (tLast < tmax) {
			x = rand() / (RAND_MAX + 1.0);
			while (x == 0) x = rand() / (RAND_MAX + 1.0);
			tLast -= log(x) / rate;
			external_exc_inputs_[i].push_back(tLast);
		}
		//cout << setprecision(15) << (double)external_exc_inputs_[i][10] << ',';
		//cout << external_exc_inputs_[i].size() << ',';
	}
	for (int i = 0; i < neuron_number_; i++) {
		rate = rates[i][1];
		tLast = external_inh_inputs_[i].front();
		while (tLast < tmax) {
			x = rand() / (RAND_MAX + 1.0);
			while (x == 0) x = rand() / (RAND_MAX + 1.0);
			tLast -= log(x) / rate;
			external_inh_inputs_[i].push_back(tLast);
		}
	}
	//cout << endl;
}

// Used in two layer network system;
// TODO: the number of sorting can be reduced;
void NeuronalNetwork::InNewSpikes(vector<vector<Spike> > & data) {
	for (int i = 0; i < neuron_number_; i++) {
		if (!data[i].empty()) {
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

void NeuronalNetwork::SetDrivingType(bool driving_type) {
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].SetDrivingType(driving_type);
	}
}

void NeuronalNetwork::UpdateNetworkState(double t, double dt) {
	if (is_con_) {
		vector<SpikeElement> T;
		double newt;
		// Creating updating pool;
		vector<int> update_list, fired_list;
		for (int i = 0; i < neuron_number_; i++) update_list.push_back(i);
		newt = SortSpikes(dym_vals_new_, update_list, fired_list, t, dt, T);
		while (newt > 0) {
			update_list.clear();
			int IND = (T.front()).index;
			fired_list.push_back(IND);
			Spike ADD_mutual;
			ADD_mutual.mode = false;
			ADD_mutual.function = (T.front()).type;
			ADD_mutual.t = t + newt + interaction_delay_;
			// erase used spiking events;
			T.erase(T.begin());
			for (int j = 0; j < neuron_number_; j++) {
				if (j == IND) {
					//neurons_[j].Fire(dym_vals_[j], t, newt);
					neurons_[j].Fire(t + newt);
				} else {
					//neurons_[j].UpdateNeuronalState(dym_vals_[j], t, newt);
					if (con_mat_[IND][j]) {
						neurons_[j].InSpike(ADD_mutual);
						update_list.push_back(j);
						// Check whether this neuron appears in the firing list T;
						for (int k = 0; k < T.size(); k ++) {
							if (j == T[k].index) {
								T.erase(T.begin() + k);
								break;
							}
						}
					}
				}
			}
			//dt -= newt;
			//t += newt;
			//T.clear();
			newt = SortSpikes(dym_vals_new_, update_list, fired_list, t, dt, T);
		}
		for (int i = 0; i < neuron_number_; i++) {
			//neurons_[i].UpdateNeuronalState(dym_vals_[i], t, dt);
			neurons_[i].UpdateNeuronalState(dym_vals_[i], dym_vals_new_[i], t + dt);
		}
	} else {
		//double spike_time;
		for (int i = 0; i < neuron_number_; i++) {
			neurons_[i].UpdateNeuronalState(dym_vals_[i], t, dt, external_exc_inputs_[i], external_inh_inputs_[i]);
			//spike_time = neurons_[i].TemporallyUpdateNeuronalState(dym_vals_[i], dym_vals_new_[i], t, dt, external_exc_inputs_[i], external_inh_inputs_[i]);
			//if (spike_time > 0) neurons_[i].Fire(t + spike_time);
			//neurons_[i].UpdateNeuronalState(dym_vals_[i], dym_vals_new_[i], t + dt);
		}
	}
}

void NeuronalNetwork::PrintCycle() {
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].GetCycle();
		cout << '\t';
	}
	cout << endl;
}
void NeuronalNetwork::OutPotential(FILEWRITE& file) {
	vector<double> potential(neuron_number_);
	for (int i = 0; i < neuron_number_; i++) {
		potential[i] = neurons_[i].GetPotential(dym_vals_[i]);
	}
	file.Write(potential);
}

void NeuronalNetwork::OutConductance(FILEWRITE& file, bool function) {
	vector<double> conductance(neuron_number_);
	if (function) {
		for (int i = 0; i < neuron_number_; i++) {
			conductance[i] = neurons_[i].GetConductance(dym_vals_[i], true);
		}
	} else {
		for (int i = 0; i < neuron_number_; i++) {
			conductance[i] = neurons_[i].GetConductance(dym_vals_[i], false);
		}
	}
	file.Write(conductance);
}

void NeuronalNetwork::OutCurrent(FILEWRITE& file) {
	vector<double> current(neuron_number_);
	for (int i = 0; i < neuron_number_; i++) {
		current[i] = neurons_[i].OutTotalCurrent(dym_vals_[i]);
	}
	file.Write(current);
}

void NeuronalNetwork::OutPartialCurrent(FILEWRITE& file, bool type) {
	vector<double> current(neuron_number_);
	for (int i = 0; i < neuron_number_; i++) {
		current[i] = neurons_[i].OutLeakyCurrent(dym_vals_[i]) + neurons_[i].OutSynapticCurrent(dym_vals_[i], type);
	}
	file.Write(current);
}

void NeuronalNetwork::SaveConMat(string connecting_matrix_file) {
	Print2D(connecting_matrix_file, con_mat_, "trunc");
}

void NeuronalNetwork::SaveNeuron(string neuron_file) {
	ofstream data;
	data.open(neuron_file.c_str(), ios::binary);
	NeuronalState add;
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].Save(add);
		data.write((char*)&add, sizeof(add));
		// data << (bool)add.type << ',';
		// data << (int)add.index << ',';
		// data << (double)add.membrane_potential_ << ',';
		// data << (double)add.ge << ',';
		// data << (double)add.gi << ',';
		// data << (double)add.remaining_refractory_time << '\n';
	}
	data.close();
}

int NeuronalNetwork::OutSpikeTrains(string path) {
	vector<vector<double> > spikes(neuron_number_);
	vector<double> add_spike_train;
	int spike_num = 0;
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].OutSpikeTrain(add_spike_train);
		spikes[i] = add_spike_train;
		spike_num += add_spike_train.size();
	}
	Print2D(path, spikes, "trunc");
	return spike_num;
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
	neurons_[i].GetConductance(dym_vals_[i], function);
}

void NeuronalNetwork::RestoreNeurons() {
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].Reset(dym_vals_[i]);
	}
	external_exc_inputs_.clear();
	external_inh_inputs_.clear();
}

void UpdateSystemState(NeuronalNetwork & pre_network, NeuronalNetwork & post_network, vector<vector<bool> > & connectivity_matrix, double t, double dt) {
	//	Update pre_network in two network system;
	pre_network.UpdateNetworkState(t, dt);
	// Get neuronal number:
	int number_pre = pre_network.GetNeuronNumber();
	int number_post = post_network.GetNeuronNumber();
	//	Transmit spiking information from pre_network to post_network;
	vector<vector<Spike> > tempPreSpikes(number_pre), tempPostSpikes(number_post);
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
	for (int i = 0; i < number_pre; i++) {
		if (!tempPreSpikes[i].empty()) {
			for (int j = 0; j < number_post; j++) {
				if (connectivity_matrix[i][j]) {
					tempPostSpikes[j].insert(tempPostSpikes[j].end(), tempPreSpikes[i].begin(), tempPreSpikes[i].end());
				}
			}
		}
	}
	post_network.InNewSpikes(tempPostSpikes);
	//	Update post_network in two network system;
	post_network.UpdateNetworkState(t, dt);
}
