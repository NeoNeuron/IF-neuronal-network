//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-02-21 16:06:04
//	Description: test program for Class NeuronalNetwork;
//***************
#include"../include/network.h"
#include <iostream>
#include <cstdio>
#include <ctime>
#include <cmath>
using namespace std;

int main() {
	clock_t start, finish;
	start = clock();
	int neuron_num = 10;
	NeuronalNetwork cells(neuron_num);
	//cells.SetS(true, 0.05);
	cells.RandNet(1, 1);
	double t = 0, dt = 0.5, tmax = 2000;
	int reps = 8;
	vector<double> add_rate = {1.5, 0};
	vector<vector<double> > PRate(neuron_num, add_rate);
	cells.InitializeNeuronalType(1, 1);
	cout << endl;
	int sampling_rate = 2;
	// prepare data file;
	string filename = "tmp/data_network_test.bin";
	size_t shape[2];
	shape[0] = reps;
	shape[1] = tmax * sampling_rate * neuron_num;
	ofstream file;
	file.open(filename.c_str(), ios::binary);
	file.write((char*)shape, 2*sizeof(size_t));
	file.close();
	// Start loop;
	for (int i = 0; i < reps; i++) {
		cells.RestoreNeurons();
		cells.SetDrivingType(true); // external type;
		cells.InitializeExternalPoissonProcess(PRate, tmax, 1);
		while (t < tmax) {
			cells.UpdateNetworkState(t, dt);
			t += dt;
			if (floor(sampling_rate * t)  == sampling_rate * t) {
				cells.OutPotential(filename);
			}
		}
		cout << i << '\n';
		t = 0;
		dt /= 2;
	}	
	int spike_num;
	spike_num = cells.OutSpikeTrains("tmp/spiketrain.csv");
	cout << "Mean firing rate: " << (double)spike_num*1000.0/tmax/neuron_num << endl;
	finish = clock();
	cout << endl;
	double ToTtime;
	ToTtime = (finish - start)*1.0 / CLOCKS_PER_SEC;
	cout << "It takes " << (double)ToTtime << " s. " << endl;	
	return 0;
}
