//******************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Description: Define class Neuron, structure Spike and NeuronState;
//	Date: 2017-02-21 16:06:30
//******************************
#include"group.h"
#include<iostream>
#include<algorithm>
#include<ctime>
#include<cmath>

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

void NeuronalNetwork::InitializeExternalPoissonProcess(bool function, double rate, double tmax, int seed) {
	vector<double> ADD;
	ADD.push_back(0);	
	srand(seed);
	double x; // temp random number;
	double tLast; // last existing data point;
	if (function == true) {
		//ofstream external_poisson;
		//external_poisson.open("external_poisson_new.txt");
		for (int i = 0; i < neuron_number_; i++) {
			external_excitatory_inputs_.push_back(ADD);
		}
		for (vector<vector<double> >::iterator it = external_excitatory_inputs_.begin(); it != external_excitatory_inputs_.end(); it++) {
			tLast = it->front();
			//external_poisson << (double)tLast << "\t";
			while (tLast < tmax) {
				x = rand() / (RAND_MAX + 1.0);
				while (log(x) < -1e12) x = rand() / (RAND_MAX + 1.0);
				tLast -= log(x) / rate;
				it->push_back(tLast);
				//external_poisson << (double)tLast << "\t";
			}
			//external_poisson << endl;
		}
		//external_poisson.close();
	} else {
		for (int i = 0; i < neuron_number_; i++) {
			external_inhibitory_inputs_.push_back(ADD);
		}
		for (vector<vector<double> >::iterator it = external_inhibitory_inputs_.begin(); it != external_inhibitory_inputs_.end(); it++) {
			tLast = it->front();
			while (tLast < tmax) {
				x = rand() / (RAND_MAX + 1.0);
				tLast -= log(x) / rate;
				it->push_back(tLast);
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

void NeuronalNetwork::LoadNetworkState(string neuron_file, string connecting_matrix_file) {
	vector<vector<double> > neuronalSetups;
	Read2DInfo(neuron_file, neuronalSetups);
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
	vector<vector<int> > connecting_matrix;
	Read2DInfo(connecting_matrix_file, connecting_matrix);
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
		mean_path = connectivity_matrix_.OutputMeanPath();
		mean_clustering_coefficient = connectivity_matrix_.OutputMeanClusteringCoefficient();
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
			for (int j = 0; j < neuron_number_; j++) {
				if (j == IND) {
					neurons_[j].Fire(t, newt - t);
				} else {
					neurons_[j].UpdateNeuronalState(t, newt - t, external_excitatory_inputs_[j], external_inhibitory_inputs_[j]);
					if (connectivity_matrix_.ReadMatrix(IND,j) != 0) {
						ADD_mutual.t = newt;
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

void NeuronalNetwork::OutPotential(vector<double>& potential) {
	potential.clear();
	potential.resize(neuron_number_);
	for (int i = 0; i < neuron_number_; i++) {
		potential[i] = neurons_[i].GetPotential();
	}
}

void NeuronalNetwork::OutTemporalParameters(vector<vector<double> > &x) {
  x.clear();
  vector<double> add(3);
  x.resize(neuron_number_, add);
  for (int i = 0; i < neuron_number_; i++) {
  	x[i][0] = neurons_[i].GetPotential();
  	x[i][1] = neurons_[i].GetConductance(true);
  	x[i][2] = neurons_[i].GetConductance(false);  
  }
}

void NeuronalNetwork::OutCurrent(vector<double> & current) {
	current.clear();
	current.resize(neuron_number_);
	for (int i = 0; i < neuron_number_; i++) {
		current[i] = neurons_[i].OutTotalCurrent();
	}
}


void NeuronalNetwork::Save(string neuron_file, string connecting_matrix_file) {
	ofstream data;
	data.open(neuron_file.c_str());
	NeuronalState add;
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].Save(add);
		data << (bool)add.type << '\t';
		data << (int)add.index << '\t';
		data << (double)add.membrane_potential_ << '\t';
		data << (double)add.ge << '\t'; 
		data << (double)add.gi << '\t';
		data << (double)add.remaining_refractory_time << endl;
	}
	data.close();
	data.open(connecting_matrix_file.c_str());
	connectivity_matrix_.OutputMatrix(data);
	data.close();
}


void NeuronalNetwork::OutSpikeTrains(vector<vector<double> > & data) {
	data.clear();
	data.resize(neuron_number_);
	vector<double> add_spike_train;
	for (int i = 0; i < neuron_number_; i++) {		
		neurons_[i].OutSpikeTrain(add_spike_train);
		data[i] = add_spike_train;
	}
}

void NeuronalNetwork::OutNewSpikes(double t, vector<vector<Spike> >& data) {
	data.clear();
	data.resize(neuron_number_);
	vector<Spike> x;
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].OutNewSpikes(t, x);
		data[i] = x;
	}
}

void NeuronalNetwork::OutNeuronType(vector<bool>& types) {
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

void NeuronalNetwork::GetMatrix(ofstream & matrix_file) {
	matrix_file.open("matFile.txt");
	connectivity_matrix_.OutputMatrix(matrix_file);
	matrix_file.close();
}

void NeuronalNetwork::RestoreNeurons() {
	for (int i = 0; i < neuron_number_; i++) {
		neurons_[i].Reset();
	}
	external_excitatory_inputs_.clear();
	external_inhibitory_inputs_.clear();
}

void Read2DInfo(string filename, vector<vector<int> >& data) {
	data.clear();
	ifstream ifile;
	const char* char_filename = filename.c_str();
	ifile.open(char_filename);
	string s;
	vector<int> add_int;
	string::size_type pos;
	while (getline(ifile, s)) {
		add_int.clear();
		pos = s.find_first_of('\t', 0);
		string ss;
		const char *sss;
		while (pos != s.npos) {
			ss = s.substr(0, pos);
			sss = ss.c_str();
			add_int.push_back(atof(sss));
			s.erase(0, pos + 1);
			ss.clear();
			pos = s.find_first_of('\t', 0);
		}
		pos = s.find_first_of('\n', 0);
		if (pos == 0) continue;
		else {
			ss = s.substr(0, pos);
			sss = ss.c_str();
			add_int.push_back(atof(sss));
		}
		data.push_back(add_int);
		s.clear();
	}
	ifile.close();
}

void Read2DInfo(string filename, vector<vector<double> >& data) {
	data.clear();
	ifstream ifile;
	const char* char_filename = filename.c_str();
	ifile.open(char_filename);
	string s;
	vector<double> add_double;
	string::size_type pos;
	while (getline(ifile, s)) {
		add_double.clear();
		pos = s.find_first_of('\t', 0);
		string ss;
		const char *sss;
		while (pos != s.npos) {
			ss = s.substr(0, pos);
			sss = ss.c_str();
			add_double.push_back(atof(sss));
			s.erase(0, pos + 1);
			ss.clear();
			pos = s.find_first_of('\t', 0);
		}
		pos = s.find_first_of('\n', 0);
		if (pos == 0) continue;
		else {
			ss = s.substr(0, pos);
			sss = ss.c_str();
			add_double.push_back(atof(sss));
		}
		data.push_back(add_double);
		s.clear();
	}
	ifile.close();
}