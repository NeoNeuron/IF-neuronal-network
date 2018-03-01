//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2018-01-30
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

//	Local field potential model [version 0.11]
//	Description: point current source model without sptial distribution;
//	DOUBLE* t_range: time period used in calculation, with unit ms, include the last point while not the first point
//	STRING current_file: total membrane current;
//	VECTOR<DOUBLE> lfp: local field potential data;
//	Return: none;
void LFP(string current_path, vector<double>& lfp, vector<int>& neuron_list, double* t_range);

#endif // _IFNET_LFP_H_
