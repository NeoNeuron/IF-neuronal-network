//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-08 21:32:55
//	Description: source file of lfp.h
//***************
#include "../include/lfp.h"
#include "../include/io.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <ctime>
using namespace std;

bool comp(const int x, const int y) {
	return x < y;
}

void Sample(vector<int> & origin_vector, vector<int> & sample_vector, int num) {
	sample_vector.clear();
	if (num == origin_vector.size() || num == 0) {
		sample_vector = origin_vector;
		return;
	} else {
		srand(time(0));
		int ind;
		sample_vector.resize(num);
		vector<int> vector_copy = origin_vector;
		for (int i = 0; i < num; i++) {
			if (vector_copy.size() == 1) ind = 0;
			else ind = rand() % (vector_copy.size() - 1);
			sample_vector[i] = vector_copy[ind];
			vector_copy.erase(vector_copy.begin() + ind, vector_copy.begin() + ind + 1);
		}
		return;
	}
}

int KeySelect(string & key, vector<neuron_type> & types, vector<int> & indices) {
	// select keys:
	indices.clear();
	if (key == "all") {
		indices.resize(types.size());
		for (int i = 0; i < types.size(); i ++) {
			indices[i] = types[i].index;
		}
	} else if (key == "exc") {
		for (vector<neuron_type>::iterator it = types.begin(); it != types.end(); it ++) {
			if (it->type == true) indices.push_back(it->index);
		}
	} else if (key == "inh") {
		for (vector<neuron_type>::iterator it = types.begin(); it != types.end(); it ++) {
			if (it->type == false) indices.push_back(it->index);
		}
	}
	return indices.size();
}

void LFP(double* t_range, vector<int> & neuron_list, string potential_path, string excitatory_conductance_path, string inhibitory_conductance_path, vector<double> &lfp) {
	// preliminary parameters;
	double sampling_rate = 32; // Unit ms: 32/ms;
	double leaky_reversal_potential = 0;
	double excitatory_reversal_potential = 14 / 3;
	double inhibitory_reversal_potential = -2 / 3;
	double leaky_conductance = 0.05;

	// Preparing time series;
	int t_begin = t_range[0] * sampling_rate; // not included
	int t_end = t_range[1] * sampling_rate; // included
	int size_of_lfp = t_end - t_begin;
	lfp.clear();
	lfp.resize(size_of_lfp);

	// Sort neuron list:
	sort(neuron_list.begin(), neuron_list.end(), comp);

	// Load potential file and conductance files;
	string s_potential, s_e_conductance, s_i_conductance;
	string ss;
	// For t = [0, t_begin];
	// cout << ">> Loading ... " << endl;
	ifstream potential_in_file, excitatory_conductance_in_file, inhibitory_conductance_in_file;
	potential_in_file.open(potential_path.c_str());
	excitatory_conductance_in_file.open(excitatory_conductance_path.c_str());
	inhibitory_conductance_in_file.open(inhibitory_conductance_path.c_str());
	for (int i = 0; i < t_begin; i ++) {
		getline(potential_in_file, s_potential);
		getline(excitatory_conductance_in_file, s_e_conductance);
		getline(inhibitory_conductance_in_file, s_i_conductance);
	}
	// For t = (t_begin, t_end]
	char cr = (char)13;
	double current_progress;
	double temp_potential, temp_e_conductance, temp_i_conductance;
	double temp_lfp;
	int neuron_list_counter;
	string::size_type pos_potential, pos_e_conductance, pos_i_conductance;
	for (int i = t_begin; i < t_end; i++) {
		getline(potential_in_file, s_potential);
		getline(excitatory_conductance_in_file, s_e_conductance);
		getline(inhibitory_conductance_in_file, s_i_conductance);
		neuron_list_counter = 0;
		temp_lfp = 0;
		for (int j = 0; j < neuron_list.back() + 1; j ++) {
			pos_potential = s_potential.find_first_of(',', 0);
			pos_e_conductance = s_e_conductance.find_first_of(',', 0);
			pos_i_conductance = s_i_conductance.find_first_of(',', 0);

			if (j == neuron_list[neuron_list_counter]) {
				ss = s_potential.substr(0, pos_potential);
				temp_potential = atof(ss.c_str());
				ss = s_e_conductance.substr(0, pos_e_conductance);
				temp_e_conductance = atof(ss.c_str());
				ss = s_i_conductance.substr(0, pos_i_conductance);
				temp_i_conductance = atof(ss.c_str());
				temp_lfp += -leaky_conductance * (temp_potential - leaky_reversal_potential) - temp_e_conductance * (temp_potential - excitatory_reversal_potential) - temp_i_conductance * (temp_potential - inhibitory_reversal_potential);
				neuron_list_counter ++;
			}
			s_potential.erase(0, pos_potential + 1);
			s_e_conductance.erase(0, pos_e_conductance + 1);
			s_i_conductance.erase(0, pos_i_conductance + 1);
			if (neuron_list_counter == neuron_list.size()) break;
		}
		lfp[i - t_begin] = temp_lfp / neuron_list.size();
		cout << cr << ">> Processing ... ";
		current_progress = (i - t_begin + 1)*100.0/size_of_lfp;
		printf("%.2f", current_progress);
		cout << "%";
	}
	potential_in_file.close();
	excitatory_conductance_in_file.close();
	inhibitory_conductance_in_file.close();
	cout << endl;
}

void LFP(double* t_range, vector<int> & neuron_list, string current_path, vector<double> &lfp) {
	// preliminary parameters;
	double sampling_rate = 32; // Unit ms: 32/ms;

	// Preparing time series;
	int t_begin = t_range[0] * sampling_rate; // not included
	int t_end = t_range[1] * sampling_rate; // included
	int size_of_lfp = t_end - t_begin;
	lfp.clear();
	lfp.resize(size_of_lfp);

	// Sort neuron list:
	sort(neuron_list.begin(), neuron_list.end(), comp);

	// Load potential file and conductance files;
	string s_current;
	string ss;
	// For t = [0, t_begin];
	// cout << ">> Loading ... " << endl;
	ifstream current_in_file;
	current_in_file.open(current_path.c_str());
	for (int i = 0; i < t_begin; i ++) {
		getline(current_in_file, s_current);
	}
	// For t = (t_begin, t_end]
	char cr = (char)13;
	double current_progress;
	double temp_current;
	double temp_lfp;
	int neuron_list_counter;
	string::size_type pos_current;
	for (int i = t_begin; i < t_end; i++) {
		getline(current_in_file, s_current);
		neuron_list_counter = 0;
		temp_lfp = 0;
		for (int j = 0; j < neuron_list.back() + 1; j ++) {
			pos_current = s_current.find_first_of(',', 0);
			if (j == neuron_list[neuron_list_counter]) {
				ss = s_current.substr(0, pos_current);
				temp_current = atof(ss.c_str());
				temp_lfp += temp_current;
				neuron_list_counter ++;
			}
			s_current.erase(0, pos_current + 1);
			if (neuron_list_counter == neuron_list.size()) break;
		}
		lfp[i - t_begin] = temp_lfp / neuron_list.size();
		cout << cr << ">> Processing ... ";
		current_progress = (i - t_begin + 1)*100.0/size_of_lfp;
		printf("%.2f", current_progress);
		cout << "%";
	}
	current_in_file.close();
	cout << endl;
}

void OutLFP(string path, vector<double>& lfp) {
	Print1D(path, "trunc", 1, lfp);
}

void OutSpikeTrain(string filename, vector<double>& spikes, double* t_range) {
	vector<double> spikes_copy = spikes;
	vector<double>::iterator it = spikes_copy.begin();
	while (it != spikes_copy.end()) {
		if (*it <= t_range[0]) it = spikes_copy.erase(it);
		else if (*it > t_range[1]) {
			spikes_copy.erase(it, spikes_copy.end());
			break;
		} else {
			*it -= t_range[0];
			it++;
		}
	}
	string path = "./data/raster/" + filename;
	Print1D(path, "trunc", 1, spikes_copy);
}
