//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2018-01-30
//	Description: source file of lfp.h
//***************
#include "../include/lfp.h"
#include "../include/io.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
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
	string buffer_potential, buffer_econd, buffer_icond;
	for (int i = t_begin; i < t_end; i++) {
		getline(potential_in_file, s_potential);
		getline(excitatory_conductance_in_file, s_e_conductance);
		getline(inhibitory_conductance_in_file, s_i_conductance);
		stringstream ipotential(s_potential);
		stringstream iecond(s_e_conductance);
		stringstream iicond(s_i_conductance);
		neuron_list_counter = 0;
		temp_lfp = 0;
		for (int j = 0; j < neuron_list.back() + 1; j ++) {
			getline(ipotential, buffer_potential, ',');
			getline(iecond, buffer_econd, ',');
			getline(iicond, buffer_icond, ',');

			if (j == neuron_list[neuron_list_counter]) {
				temp_potential = atof(buffer_potential.c_str());
				temp_e_conductance = atof(buffer_econd.c_str());
				temp_i_conductance = atof(buffer_icond.c_str());
				temp_lfp += -leaky_conductance * (temp_potential - leaky_reversal_potential) - temp_e_conductance * (temp_potential - excitatory_reversal_potential) - temp_i_conductance * (temp_potential - inhibitory_reversal_potential);
				neuron_list_counter ++;
			}
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
	// For t = [0, t_begin];
	// cout << ">> Loading ... " << endl;
	ifstream current_in_file;
	current_in_file.open(current_path.c_str());
	for (int i = 0; i < t_begin; i ++) {
		getline(current_in_file, s_current);
	}
	// For t = (t_begin, t_end]
	char cr = (char)13;
	int current_progress = 0;
	double temp_lfp;
	size_t neuron_list_counter;
	string buffer;
	for (int i = t_begin; i < t_end; i++) {
		getline(current_in_file, s_current);
		stringstream icurrent(s_current);
		neuron_list_counter = 0;
		temp_lfp = 0;
		for (int j = 0; j < neuron_list.back() + 1; j ++) {
			getline(icurrent, buffer, ',');
			if (j == neuron_list[neuron_list_counter]) {
				temp_lfp += atof(buffer.c_str());
				neuron_list_counter ++;
			}
			if (neuron_list_counter == neuron_list.size()) break;
		}
		lfp[i - t_begin] = temp_lfp / neuron_list.size();
		if (floor((i - t_begin + 1)*100.0/size_of_lfp - current_progress) >= 1) {
			current_progress ++;
			cout << cr << ">> Processing ... ";
			printf("%d", current_progress);
			cout << "%";
		}
	}
	current_in_file.close();
	cout << endl;
}
