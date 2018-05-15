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

void InitializeBinFile(string filename, size_t* shape) {
	ofstream file;
	file.open(filename.c_str(), ios::binary);
	file.write((char*)shape, 2*sizeof(size_t));
	file.close();
}

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
  cout << ">> [Config.ini]:\n#####\n";
	PrintConfig(m_map_config);
	cout << "#####\n";
	// load neuron number;
	int neuron_number = atoi(m_map_config["NeuronNumber"].c_str());
	NeuronalNetwork net(neuron_number);
	// load connecting mode;
	int connecting_mode = atoi(m_map_config["ConnectingMode"].c_str());
	if (connecting_mode == 0) { // External connectivity matrix;
		net.LoadConnectivityMat("./doc/externalMat.csv");
	} else if (connecting_mode == 1) {
		int connecting_density = atoi(m_map_config["ConnectingDensity"].c_str());
		net.SetConnectingDensity(connecting_density);
		double rewiring_probability = atof(m_map_config["RewiringProbability"].c_str());
		int rewiring_seed = atoi(m_map_config["NetSeed"].c_str());
		// Generate networks;
		net.Rewire(rewiring_probability, rewiring_seed, true);
	} else if (connecting_mode == 2) {
		net.RandNet(atof(m_map_config["RandomProbability"].c_str()), atoi(m_map_config["NetSeed"].c_str()));
	}
	// Set interneuronal coupling strength;
	net.SetS(true, atof(m_map_config["SynapticStrengthExcitatory"].c_str()));
	net.SetDelay(0.0);
	// initialize the network;
	net.InitializeNeuronalType(atof(m_map_config["TypeProbability"].c_str()), atoi(m_map_config["TypeSeed"].c_str()));
	cout << "in the network." << endl;

	double maximum_time = atof(m_map_config["MaximumTime"].c_str());
	bool driving_type, fwd_type;
	istringstream(m_map_config["DrivingType"]) >> boolalpha >> driving_type;
	istringstream(m_map_config["HomoDriving"]) >> boolalpha >> fwd_type;
	vector<double> fwd_rates(neuron_number);
	net.SetDrivingType(driving_type);
	if (fwd_type) {
		double rate_exc = atof(m_map_config["DrivingRateExcitatory"].c_str());
		double rate_inh = atof(m_map_config["DrivingRateInhibitory"].c_str());
		vector<bool> neuron_types;
		net.GetNeuronType(neuron_types);
		for (int i = 0; i < neuron_number; i ++) {
			if (neuron_types[i]) {
				fwd_rates[i] = rate_exc;
			} else {
				fwd_rates[i] = rate_inh;
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
		Read1D("./doc/fwd_setting.csv", fwd_rates, 0, 1);
		if (driving_type) {
			net.InitializeExternalPoissonProcess(fwd_rates, maximum_time, atoi(m_map_config["ExternalDrivingSeed"].c_str()));
		} else {
			srand(time(NULL));
			net.InitializeInternalPoissonRate(fwd_rates);
		}
	}

	// Set driving strength;
	net.SetF(atof(m_map_config["DrivingStrength"].c_str()));

	// SETUP DYNAMICS:
	double t = 0, dt = atof(m_map_config["TimingStep"].c_str()), tmax = maximum_time;
	double recording_rate = 2;
	// Define the shape of data;
	size_t shape[2];
	shape[0] = tmax * recording_rate;
	shape[1] = neuron_number;
	// Define file path for output data;
	string V_path = dir + "V.bin";
	string I_path = dir + "I.bin";
	//string GE_path = dir + "GE.bin";
	// Initialize files:
	InitializeBinFile(V_path, shape);
	InitializeBinFile(I_path, shape);
	//InitializeBinFile(GE_path, shape);

	// double progress;
	while (t < tmax) {
		net.UpdateNetworkState(t, dt);
		t += dt;
		// Output temporal data;
		if (abs(recording_rate*t - floor(recording_rate*t)) == 0) {
			net.OutPotential(V_path);
			net.OutCurrent(I_path);
		}
		//net.OutConductance(GE_path, true);
	}

	string neuron_path, mat_path;
	neuron_path = dir + "neuron.bin";
	mat_path = dir + "mat.csv";
	net.Save(neuron_path, mat_path);

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
