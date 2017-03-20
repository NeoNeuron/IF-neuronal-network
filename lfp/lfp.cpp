//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-08 21:32:55
//	Description: source file of lfp.h
//***************
#include "lfp.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <ctime>

using namespace std;

void Sample(vector<int> & origin_vector, vector<int> & sample_vector, int num) {
	if (origin_vector.size() == num) return;
	else {
		srand(time(0));
		int ind;
		sample_vector.clear();
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


void ReadLine(string file_name, int line_index, vector<int> &output) {
	// open data file;
	const char* char_file_name = file_name.c_str();
	ifstream ifile;
	ifile.open(char_file_name);
	// prepare input file stream;
	string s;
	output.clear();
	int getline_counter = 0;
	while (getline(ifile, s)) {
		if (getline_counter == line_index) {
			string::size_type pos = s.find_first_of('\t', 0);
			string ss;
			while (pos != s.npos) {
				ss = s.substr(0, pos);
				const char *sss;
				sss = ss.c_str();
				output.push_back(atof(sss));
				s.erase(0, pos + 1);
				ss.clear();
				pos = s.find_first_of('\t', 0);
			}
			pos = s.find_first_of('\n', 0);
			if (pos != 0) {
				ss = s.substr(0, pos);
				const char *sss;
				sss = ss.c_str();
				output.push_back(atof(sss));
			}
			break;
		}
		getline_counter ++;
	}
	ifile.close();
}

void ReadLine(string file_name, int line_index, vector<double> &output) {
	// open data file;
	const char* char_file_name = file_name.c_str();
	ifstream ifile;
	ifile.open(char_file_name);
	// prepare input file stream;
	string s;
	output.clear();
	int getline_counter = 0;
	while (getline(ifile, s)) {
		if (getline_counter == line_index) {
			string::size_type pos = s.find_first_of('\t', 0);
			string ss;
			while (pos != s.npos) {
				ss = s.substr(0, pos);
				const char *sss;
				sss = ss.c_str();
				output.push_back(atof(sss));
				s.erase(0, pos + 1);
				ss.clear();
				pos = s.find_first_of('\t', 0);
			}
			pos = s.find_first_of('\n', 0);
			if (pos != 0) {
				ss = s.substr(0, pos);
				const char *sss;
				sss = ss.c_str();
				output.push_back(atof(sss));
			}
			break;
		}
		getline_counter ++;
	}
	ifile.close();
}

void ReadColumn(string file_name, int column_index, int num_column, vector<int> &output) {
	// open data file;
	const char* char_file_name = file_name.c_str();
	ifstream ifile;
	ifile.open(char_file_name);
	// prepare input file stream;
	string s;
	string ss;
	string::size_type pos;
	int column_counter;
	output.clear();
	while (getline(ifile, s)) {
		for (int i = 0; i < num_column; i++) {
			column_counter = 0;
			pos = s.find_first_of('\t', 0);
			if (column_counter == column_index) {
				ss = s.substr(0, pos);
				const char *sss;
				sss = ss.c_str();
				output.push_back(atof(sss));
				break;
			}
			s.erase(0, pos + 1);
			column_counter ++;
		}
	}
	ifile.close();
}

void ReadLines(string file_name, vector<int> &line_index, vector<vector<double> > &data) {
	// open data file;
	const char* char_file_name = file_name.c_str();
	ifstream ifile;
	ifile.open(char_file_name);
	// prepare input file stream;
	string s;
	vector<double> add_double;
	data.clear();
	int getline_counter = -1;
	int line_index_i = 0;
	int line_index_max = line_index.size();
	while (getline(ifile, s)) {
		getline_counter ++;
		//cout << getline_counter << endl;
		if (line_index[line_index_i] == getline_counter) {
			add_double.clear();
			string::size_type pos = s.find_first_of('\t', 0);
			string ss;			
			//int counter = 0;
			while (pos != s.npos) {
				ss = s.substr(0, pos);
				const char *sss;
				sss = ss.c_str();
				add_double.push_back(atof(sss));
				s.erase(0, pos + 1);
				ss.clear();
				pos = s.find_first_of('\t', 0);
				//counter ++;
				//cout << counter << endl;
			}
			pos = s.find_first_of('\n', 0);
			if (pos != 0) {
				ss = s.substr(0, pos);
				const char *sss;
				sss = ss.c_str();
				add_double.push_back(atof(sss));
			}
			data.push_back(add_double);
			s.clear();
			line_index_i ++;
		}
		if (line_index_i == line_index_max) break;
	}
	ifile.close();
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


void LFP(double* t_range, vector<int> & neuron_list, string potential_filename, string excitatory_conductance_filename, string inhibitory_conductance_filename, vector<double> &lfp) {
	// preliminary parameters;
	double sampling_rate = 32; // Unit ms: 32/ms;
	double leaky_reversal_potential = 0;
	double excitatory_reversal_potential = 14 / 3;
	double inhibitory_reversal_potential = -2 / 3;
	double leaky_conductance = 0.05;
	int total_neuron_number = 100;

	// Preparing time series;
	int t_begin = t_range[0] * sampling_rate; // not included
	int t_end = t_range[1] * sampling_rate; // included
	int size_of_lfp = t_end - t_begin;
	lfp.clear();
	lfp.resize(size_of_lfp);

	// Load potential file and conductance files;
	const char* char_potential_filename = potential_filename.c_str();
	const char* char_excitatory_conductance_filename = excitatory_conductance_filename.c_str();
	const char* char_inhibitory_conductance_filename =inhibitory_conductance_filename.c_str();
	string s_potential, s_e_conductance, s_i_conductance;
	string ss;
	// For t = [0, t_begin];
	cout << ">> Loading ... " << endl;
	ifstream potential_in_file, excitatory_conductance_in_file, inhibitory_conductance_in_file;
	potential_in_file.open(char_potential_filename);
	excitatory_conductance_in_file.open(char_excitatory_conductance_filename);
	inhibitory_conductance_in_file.open(char_inhibitory_conductance_filename);
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
		for (int j = 0; j < total_neuron_number; j ++) {
			pos_potential = s_potential.find_first_of('\t', 0);
			pos_e_conductance = s_e_conductance.find_first_of('\t', 0);
			pos_i_conductance = s_i_conductance.find_first_of('\t', 0);
			
			if (j == neuron_list[neuron_list_counter]) {	
				ss = s_potential.substr(0, pos_potential);
				const char* sss1 = ss.c_str();
				temp_potential = atof(sss1);
				ss = s_e_conductance.substr(0, pos_e_conductance);
				const char* sss2 = ss.c_str();
				temp_e_conductance = atof(sss2);

				ss = s_i_conductance.substr(0, pos_i_conductance);
				const char* sss3 = ss.c_str();
				temp_i_conductance = atof(sss3);
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

void OutputLFP(vector<double> &lfp, string filename) {
	ofstream output;
	 const char* char_filename = filename.c_str();
	output.open(char_filename);
	for (vector<double>::iterator it = lfp.begin(); it != lfp.end(); it++) {
		output << setprecision(20) << (double)*it;
		if (it != lfp.end() - 1) output << endl;
	}
	output.close();
}

void OutputSpikeTrain(double* t_range, vector<double> &spikes, string filename) {
	ofstream output;
	 const char* char_filename = filename.c_str();
	output.open(char_filename);
	for (vector<double>::iterator it = spikes.begin(); it != spikes.end(); it++) {
		if (*it > t_range[0] and *it <= t_range[1]) {
			output << (double)*it - t_range[0] << endl;
		}
		if (*it > t_range[1]) break;
	}
	output.close();
}
