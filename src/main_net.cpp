//*************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-13 15:07:52
//	Description: test program for multi-network simulation;
//*************************

#include "../include/network.h"
#include "../include/get-config.h"
#include "../include/io.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <sstream>
#include <stdexcept>
using namespace std;

//	Simulation program for single network system;
//	arguments:
//	argv[1] = Outputing directory for neur al data;
int main(int argc, const char* argv[]) {
	if (argc != 2) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// 	Setup directory for output files;
	//	it must be existing dir;
	string dir;
	dir = argv[1];

	// Loading config.ini:
	string net_config_path = "./doc/config_net.ini";
  map<string, string> m_map_config;
  ReadConfig(net_config_path,m_map_config);
  cout << ">> [Config.ini]:\n#####";
	PrintConfig(m_map_config);
	cout << "\n#####\n";
	// load neuron number;
	int neuron_number = atoi(m_map_config["NeuronNumber"].c_str());
	NeuronalNetwork net(neuron_number);
	// load connecting mode;
	int connecting_mode = atoi(m_map_config["ConnectingMode"].c_str());
	if (connecting_mode == 0) { // External connectivity matrix;
		net.LoadConnectivityMat("./data/sampleMat.csv");
	} else {
		int connecting_density = atoi(m_map_config["ConnectingDensity"].c_str());
		net.SetConnectingDensity(connecting_density);
		int rewiring_probability = atof(m_map_config["RewiringProbability"].c_str());
		int rewiring_seed = atoi(m_map_config["RewiringSeed"].c_str());
		// Generate networks;
		net.Rewire(rewiring_probability, rewiring_seed, true);
	}
	// Set interneuronal coupling strength;
	net.SetS(true, atof(m_map_config["SynapticStrengthExcitatory"].c_str()));
	// initialize the network;
	net.InitializeNeuronalType(atof(m_map_config["TypeProbability"].c_str()), atoi(m_map_config["TypeSeed"].c_str()));
	cout << "in the network." << endl;

	double maximum_time = atof(m_map_config["MaximumTime"].c_str());
	bool driving_type, fwd_type;
	istringstream(m_map_config["DrivingType"]) >> boolalpha >> driving_type;
	istringstream(m_map_config["HomoDriving"]) >> boolalpha >> fwd_type;
	vector<vector<double> > fwd_rates;
	net.SetDrivingType(driving_type);
	if (fwd_type) {
		double rate_exc = atof(m_map_config["DrivingRateExcitatory"].c_str());
		double rate_inh = atof(m_map_config["DrivingRateInhibitory"].c_str());
		vector<bool> neuron_types;
		net.GetNeuronType(neuron_types);
		for (int i = 0; i < neuron_number; i ++) {
			if (neuron_types[i]) {
				fwd_rates[i].push_back(rate_exc);
				fwd_rates[i].push_back(0.0);
			} else {
				fwd_rates[i].push_back(rate_inh);
				fwd_rates[i].push_back(0.0);
			}
		}
		if (driving_type) {
			net.InitializeExternalPoissonProcess(fwd_rates, maximum_time, atoi(m_map_config["ExternalDrivingSeed"].c_str()));
		} else {
			srand(time(NULL));
			net.InitializeInternalPoissonRate(fwd_rates);
		}
	} else {
		// import the data file of feedforward driving rate:
		Read2D("./doc/fwd_setting.csv", fwd_rates);
		if (driving_type) {
			net.InitializeExternalPoissonProcess(fwd_rates, maximum_time, atoi(m_map_config["ExternalDrivingSeed"].c_str()));
		} else {
			srand(time(NULL));
			net.InitializeInternalPoissonRate(fwd_rates);
		}
	}

	// Set driving strength;
	net.SetF(true, atof(m_map_config["DrivingStrength"].c_str()));

	// SETUP DYNAMICS:
	double t = 0, dt = atof(m_map_config["TimingStep"].c_str()), tmax = maximum_time;
	// Define file path for output data;
	// string V_path = dir + "V.csv";
	string I_path = dir + "I.csv";
	// Initialize files:
	ofstream V, I;
	// V.open(V_path.c_str());
	// V.close();

	I.open(I_path.c_str());
	I.close();

	// char cr = (char)13;
	// double progress;
	while (t < tmax) {
		net.UpdateNetworkState(t, dt);
		t += dt;
		// Output temporal data;
		// net.OutPotential(V_path);
		net.OutCurrent(I_path);

		// progress = t * 100.0 / tmax;
		// cout << cr;
		// printf(">> Processing ... %6.2f", progress);
		// cout << "%";
	}
	// cout << endl;

	// string neuron_path, mat_path;
	// neuron_path = dir + "neuron.csv";
	// mat_path = dir + "mat.csv";
	// net.Save(neuron_path, mat_path);

	// OUTPUTS:
	string raster_path = dir + "raster.csv";
	net.OutSpikeTrains(raster_path);

	finish = clock();

	// COUNTS:
	double ToTtime;
	ToTtime = (finish - start)*1.0 / CLOCKS_PER_SEC;
	printf(">> It takes %.2fs\n", ToTtime);
	return 0;
}
