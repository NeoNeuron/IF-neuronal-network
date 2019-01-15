//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-02-21 16:06:04
//	Description: test program for Class NeuronalNetwork;
//***************
#include"../include/network.h"
#include"../include/get-config.h"
#include <iostream>
#include <cstdio>
#include <ctime>
#include <cmath>
using namespace std;

int main() {
	clock_t start, finish;
	start = clock();
	// Load config file
  map<string, string> m_map_config;
  ReadConfig("./test/config_net.ini", m_map_config);
  cout << ">> [Config.ini]:\n#####\n";
	PrintConfig(m_map_config);
	cout << "#####\n";
	int neuron_num = atoi(m_map_config["NeuronNumber"].c_str());
	NeuronalNetwork cells(m_map_config["NeuronType"], neuron_num);
	cells.InitializeNeuronalType(m_map_config);
	cells.InitializeConnectivity(m_map_config);
	cells.InitializeSynapticStrength(m_map_config);
	cells.InitializeSynapticDelay(m_map_config);
	double t = 0;
	double dt = atof(m_map_config["StartingTimingStep"].c_str());
	double tmax = atof(m_map_config["MaximumTime"].c_str());
	int reps = atoi(m_map_config["Reps"].c_str());
	double sampling_rate = 1.0 / atof(m_map_config["SamplingTimingStep"].c_str());
	// prepare data file;
	FILEWRITE file("./tmp/data_network_test.bin", "trunc");
	size_t shape[2];
	shape[0] = reps;
	shape[1] = tmax * sampling_rate * neuron_num;
	file.SetSize(shape);
	int spike_num;
	vector<vector<double> > spike_trains;
	vector<double> add_spike_train;
	// Start loop;
	for (int i = 0; i < reps; i++) {
		cells.RestoreNeurons();
		cells.InitializePoissonGenerator(m_map_config);
		while (t < tmax) {
			cells.UpdateNetworkState(t, dt);
			t += dt;
			if (floor(sampling_rate * t) == sampling_rate * t) {
				cells.OutPotential(file);
			}
		}
		spike_num = cells.OutSpikeTrains(spike_trains);
		add_spike_train.clear();
		for (int i = 0; i < spike_trains.size(); i ++) {
			add_spike_train.insert(add_spike_train.end(), spike_trains[i].begin(), spike_trains[i].end());
		}
		Print1D("./tmp/data_network_raster.csv", add_spike_train, "app", 0);
		printf("[-] dt = %.2e s\tmean firing rate = %.2f Hz\n", dt, spike_num*1000.0/tmax/neuron_num);
		t = 0;
		dt /= 2;
	}	
	finish = clock();
	cout << endl;
	double ToTtime;
	ToTtime = (finish - start)*1.0 / CLOCKS_PER_SEC;
	cout << "It takes " << (double)ToTtime << " s. " << endl;	
	return 0;
}
