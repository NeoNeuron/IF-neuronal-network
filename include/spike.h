#ifndef _IFNET_SPIKE_H_
#define _IFNET_SPIKE_H_
#include <string>
#include <vector>
using namespace std;

//	Convert double spike train to binary sequence;
//	VECTOR<DOUBLE> spikes: original spike train;
//	VECTOR<BOOL> binary_spikes: binary sequence; true for having spike; false for no spike;
//	DOUBLE tmax: maximum time of time range;
//	DOUBLE dt: the size of time step that a binary value considered;
//	Return: none;
void Spike2Bool(vector<double> & spikes, vector<bool> & binary_spikes, double tmax, double dt);

// Truncate spikes withn selected time range;
void Truncate(vector<double>& spikes, double* t_range);

// 	Output spike train;
//	STRING filename: output filename;
//	VECTOR<DOUBLE> spikes: the total spike train of target neuron;
//	Return: none;
void OutSpikeTrain(string path, vector<bool>& spikes);

#endif
