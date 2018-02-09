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

void LFP(string current_path, vector<double>& lfp, vector<int>& neuron_list, double* t_range) {
	// preliminary parameters;
	double sampling_rate = 32; // Unit ms: 32/ms;

	// Preparing time series;
	int t_begin = t_range[0] * sampling_rate; // not included
	int t_end = t_range[1] * sampling_rate; // included
	int size_of_lfp = t_end - t_begin;
	lfp.clear();
	lfp.resize(size_of_lfp);
	// Load current file;
	ifstream current_in_file;
	current_in_file.open(current_path.c_str(), ios::binary);
	size_t shape[2];
	current_in_file.read((char*)&shape, 2*sizeof(size_t));

	vector<int> diff_list(neuron_list.size() + 1);
	// Sort neuron list:
	if (neuron_list.size() == 1) {
		diff_list[0] = neuron_list[0];
		diff_list[1] = shape[1] - neuron_list[0] - 1;
	} else {
		sort(neuron_list.begin(), neuron_list.end(), comp);
		for (int i = 0; i < diff_list.size(); i ++) {
			if (i == 0) diff_list[i] = neuron_list[0];
			else if (i == diff_list.size() - 1) diff_list[i] = shape[1] - neuron_list[i - 1] - 1;
			else diff_list[i] = neuron_list[i] - neuron_list[i - 1] - 1;
		}
	}
	// For t = [0, t_begin];
	current_in_file.seekg(shape[1]*t_begin*sizeof(double), current_in_file.cur);
	// For t = (t_begin, t_end]
	double temp_lfp;
	double buffer;
	for (int i = t_begin; i < t_end; i++) {
		temp_lfp = 0;
		for (int j = 0; j < diff_list.size() - 1; j ++) {
			current_in_file.seekg(diff_list[j]*sizeof(double), current_in_file.cur);
			current_in_file.read((char*)&buffer, sizeof(double));
			temp_lfp += buffer;
		}
		current_in_file.seekg(*(diff_list.end() - 1)*sizeof(double), current_in_file.cur);
		lfp[i - t_begin] = temp_lfp / neuron_list.size();
	}
	current_in_file.close();
}
