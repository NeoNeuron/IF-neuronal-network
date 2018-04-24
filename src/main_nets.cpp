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

//	Simulation program for two network system;
//	arguments:
//	argv[1] = Outputing directory for neural data;
int main(int argc, const char* argv[]) {
	if (argc != 2) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// 	Setup directory for output files;
	//	it must be existing dir;
	string dir;
	dir = argv[1];

	// Loading config.ini:
	string net_config_path = "./doc/config_nets.ini";
  map<string, string> m_map_config;
  ReadConfig(net_config_path,m_map_config);
  cout << ">> [Configs]:" << endl;
	PrintConfig(m_map_config);
	cout << endl;

	int preNetNum, postNetNum, preNetDensity, postNetDensity;
	double pre_rewiring_probability, post_rewiring_probability;
	int pre_rewiring_seed, post_rewiring_seed;
	// neuron number:
	preNetNum = atoi(m_map_config["PreNetNeuronNumber"].c_str());
	postNetNum = atoi(m_map_config["PostNetNeuronNumber"].c_str());
	// connecting density:
	preNetDensity = atoi(m_map_config["PreNetConnectingDensity"].c_str());
	postNetDensity = atoi(m_map_config["PostNetConnectingDensity"].c_str());
	// rewiring probability and rewiring seed;
	pre_rewiring_probability = atof(m_map_config["PreNetRewiringProbability"].c_str());
	post_rewiring_probability = atof(m_map_config["PostNetRewiringProbability"].c_str());
	pre_rewiring_seed = atoi(m_map_config["PreNetRewiringSeed"].c_str());
	post_rewiring_seed = atoi(m_map_config["PostNetRewiringSeed"].c_str());
	// Initialize networks;
	NeuronalNetwork preNet(preNetNum), postNet(postNetNum);
	preNet.SetConnectingDensity(preNetDensity);
	postNet.SetConnectingDensity(postNetDensity);
	preNet.Rewire(pre_rewiring_probability, pre_rewiring_seed, true);
	postNet.Rewire(post_rewiring_probability, post_rewiring_seed, true);
	double maximum_time = atof(m_map_config["MaximumTime"].c_str());
	// Initialize neuronal types;
	preNet.InitializeNeuronalType(atof(m_map_config["PreNetTypeProbability"].c_str()), atoi(m_map_config["PreNetTypeSeed"].c_str()));
	cout << "in pre-network." << endl;
	postNet.InitializeNeuronalType(atof(m_map_config["PostNetTypeProbability"].c_str()), atoi(m_map_config["PostNetTypeSeed"].c_str()));
	cout << "in post-network." << endl;
	// Initialize driving type;
	bool preType, postType;
	istringstream(m_map_config["PreNetDrivingType"]) >> boolalpha >> preType;
	istringstream(m_map_config["PostNetDrivingType"]) >> boolalpha >> postType;
	preNet.SetDrivingType(preType);
	postNet.SetDrivingType(postType);
	// Initialize external driving;
	double prenet_rate_exc = atof(m_map_config["PreNetDrivingRateExcitatory"].c_str());
	double prenet_rate_inh = atof(m_map_config["PreNetDrivingRateInhibitory"].c_str());
	double postnet_rate_exc = atof(m_map_config["PostNetDrivingRateExcitatory"].c_str());
	double postnet_rate_inh = atof(m_map_config["PostNetDrivingRateInhibitory"].c_str());
	vector<bool> prenet_types, postnet_types;
	vector<double> prenet_fwd_rates(preNetNum), postnet_fwd_rates(postNetNum);
	preNet.GetNeuronType(prenet_types);
	postNet.GetNeuronType(postnet_types);
	for (int i = 0; i < preNetNum; i ++) {
		if (prenet_types[i]) {
			prenet_fwd_rates[i] = prenet_rate_exc;
		} else {
			prenet_fwd_rates[i] = prenet_rate_inh;
		}
	}
	for (int i = 0; i < postNetNum; i ++) {
		if (postnet_types[i]) {
			postnet_fwd_rates[i] = postnet_rate_exc;
		} else {
			postnet_fwd_rates[i] = postnet_rate_inh;
		}
	}
	if (preType) {
		preNet.InitializeExternalPoissonProcess(prenet_fwd_rates, maximum_time, atoi(m_map_config["PreNetExternalDrivingSeed"].c_str()));
	} else preNet.InitializeInternalPoissonRate(prenet_fwd_rates);
	if (postType) {
		postNet.InitializeExternalPoissonProcess(postnet_fwd_rates, maximum_time, atoi(m_map_config["PostNetExternalDrivingSeed"].c_str()));
	} else postNet.InitializeInternalPoissonRate(postnet_fwd_rates);
	// Set feedforward driving strength;
	preNet.SetF(atof(m_map_config["PreNetDrivingStrength"].c_str()));
	postNet.SetF(atof(m_map_config["PostNetDrivingStrength"].c_str()));

	// Set synaptic strength of presynaptic neural network;
	preNet.SetS(true, atof(m_map_config["PreNetSynapticStrengthExcitatory"].c_str()));
	postNet.SetS(true, atof(m_map_config["PostNetSynapticStrengthExcitatory"].c_str())/postNetDensity);

	// Setup connectivity matrix;
	vector<vector<bool> > conMat(preNetNum, vector<bool>(postNetNum));
	double conP = atof(m_map_config["ConnectingProbability"].c_str());
	srand(atoi(m_map_config["ConnectingSeed"].c_str()));
	for (int i = 0; i < preNetNum; i++) {
		for (int j = 0; j < postNetNum; j++) {
			if (rand() / (RAND_MAX*1.0) < conP) conMat[i][j] = true;
			else conMat[i][j] = false;
		}
	}
	// Output connectivity matrix;
	string conMat_path = dir + "conMat.txt";
	Print2D(conMat_path, conMat, "trunc");

	// SETUP DYNAMICS:
	double t = 0, dt = atof(m_map_config["TimingStep"].c_str()), tmax = maximum_time;
	// Define the shape of data;
	size_t preNetshape[2], postNetshape[2];
	preNetshape[0] = tmax / dt;
	postNetshape[0] = tmax /dt;
	preNetshape[1] = preNetNum;
	postNetshape[1] = postNetNum;

	// Define file path for output data;
	string preV = dir + "preV.bin";
	string postV = dir + "postV.bin";
	string preI = dir + "preI.bin";
	string postI = dir + "postI.bin";
	// Initialize files:
	InitializeBinFile(preV, preNetshape);
	InitializeBinFile(postV, postNetshape);
	InitializeBinFile(preI, preNetshape);
	InitializeBinFile(postI, postNetshape);

	// char cr = (char)13;
	int progress;
	while (t < tmax) {
		UpdateSystemState(preNet, postNet, conMat, t, dt);
		t += dt;
		// Output temporal data;
		preNet.OutPotential(preV);
		preNet.OutCurrent(preI);
		postNet.OutPotential(postV);
		postNet.OutCurrent(postI);

		if (floor(t * 10 / tmax) > progress) {
			progress = floor(t * 10/tmax);
			printf(">> Processing ... %d0%\n", progress);
		}
		// cout << cr;
	}
	// cout << endl;

	// string preNeuron_path, preMat_path;
	// preNeuron_path = dir + "preNeuron.csv";
	// preMat_path = dir + "preMat.csv";
	// preNet.Save(preNeuron_path, preMat_path);
	// string postNeuron_path, postMat_path;
	// postNeuron_path = dir + "postNeuron.csv";
	// postMat_path = dir + "postMat.csv";
	// postNet.Save(postNeuron_path, postMat_path);

	// OUTPUTS:
	string pre_raster_path = dir + "rasterPre.csv";
	string post_raster_path = dir + "rasterPost.csv";
	preNet.OutSpikeTrains(pre_raster_path);
	postNet.OutSpikeTrains(post_raster_path);

	finish = clock();

	// COUNTS:
	printf(">> It takes %.2fs\n", (finish - start)*1.0 / CLOCKS_PER_SEC);
	return 0;
}
