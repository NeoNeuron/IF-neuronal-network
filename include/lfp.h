//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-08 21:32:46
//	Description: Define model of local field potential(LFP);
//***************
#ifndef _IFNET_LFP_H_
#define _IFNET_LFP_H_

#include <string>
#include <vector>

using namespace std;

struct neuron_type{
	int index;
	bool type;
};

void Sample(vector<int> & origin_vector, vector<int> &sample_vector, int num);

//	Key Selection:
//	Choose the specific portion of neurons in post-network, according to key;
//	Return number of selected neurons;
int KeySelect(string & key, vector<neuron_type> & type, vector<int> & indices);

//	Local field potential model [version 0.10]
//	Description: point current source model without sptial distribution;
//	DOUBLE* t_range: time period used in calculation, with unit ms, include the last point while not the first point
//	STRING potential_file: membrane potential;
//	STRING excitatory_conductance_file;
//	STRING inhibitory_conductance_file;
//	VECTOR<DOUBLE> lfp: local field potential data;
//	Return: none;
void LFP(double* t_range, int total_neuron_number, vector<int> & neuron_list, string potential_path, string excitatory_conductance_path, string inhibitory_conductance_path, vector<double> &lfp);

//	Local field potential model [version 0.11]
//	Description: point current source model without sptial distribution;
//	DOUBLE* t_range: time period used in calculation, with unit ms, include the last point while not the first point
//	STRING current_file: total membrane current;
//	VECTOR<DOUBLE> lfp: local field potential data;
//	Return: none;
void LFP(double* t_range, int total_neuron_number, vector<int> & neuron_list, string current_path, vector<double> &lfp);

// 	Output LFP;
//	Description: output LFP data to a given file;
//	STRING path: output file name;
//	Return: none;
void OutLFP(string path, vector<double>& lfp);

// 	Output spike train;
//	Description: output the spike train of chosen neuron within chosen time range; spiking time points are rearranged which starts from 0;
//	DOUBLE* t_range: time period used in calculation, with unit ms;
//	STRING path: output file name;
//	Return: none;
void OutSpikeTrain(string path, vector<double>& spikes, double* t_range);

#endif // _IFNET_LFP_H_
