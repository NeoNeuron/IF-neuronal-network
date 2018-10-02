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
	NeuronalNetwork cells(neuron_num);
	cells.InitializeConnectivity(m_map_config);
	cells.InitializeSynapticStrength(m_map_config);
	double t = 0;
	double dt = atof(m_map_config["StartingTimingStep"].c_str());
	double tmax = atof(m_map_config["MaximumTime"].c_str());
	int reps = atoi(m_map_config["Reps"].c_str());
	double pre = atof(m_map_config["DrivingRateE"].c_str());
	double pri = atof(m_map_config["DrivingRateI"].c_str());
	double pse = atof(m_map_config["DrivingStrengthE"].c_str());
	double psi = atof(m_map_config["DrivingStrengthI"].c_str());
	vector<vector<double> > Pdriving(neuron_num, vector<double>{pre, pri, pse, psi});
	cells.InitializeNeuronalType(atof(m_map_config["TypeProbability"].c_str()), atoi(m_map_config["TypeSeed"].c_str()));
	cout << endl;
	double sampling_rate = 1.0 / atof(m_map_config["SamplingTimingStep"].c_str());
	// prepare data file;
	FILEWRITE file("./tmp/data_network_test.bin", "trunc");
	size_t shape[2];
	shape[0] = reps;
	shape[1] = tmax * sampling_rate * neuron_num;
	file.SetSize(shape);
	int spike_num;
	// Start loop;
	for (int i = 0; i < reps; i++) {
		cells.RestoreNeurons();
		cells.SetDrivingType(true); // external type;
		cells.InitializeExternalPoissonProcess(Pdriving, tmax, 3);
		while (t < tmax) {
			cells.UpdateNetworkState(t, dt);
			t += dt;
			if (floor(sampling_rate * t)  == sampling_rate * t) {
				cells.OutPotential(file);
			}
		}
		spike_num = cells.OutSpikeTrains("tmp/spiketrain.csv");
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
