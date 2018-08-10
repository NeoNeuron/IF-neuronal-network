//*************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2018-07-29
//	Description: test program for multi-network simulation;
//*************************

#include "../include/network.h"
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

	int preNetNum, postNetNum;
	// neuron number:
	preNetNum = atoi(m_map_config["PreNetNeuronNumber"].c_str());
	postNetNum = atoi(m_map_config["PostNetNeuronNumber"].c_str());
	// Initialize networks;
	NeuronalNetwork preNet(preNetNum), postNet(postNetNum);
	preNet.InitializeConnectivity(m_map_config, "PreNet");
	postNet.InitializeConnectivity(m_map_config, "PostNet");

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
	vector<vector<double> > prenet_fwd_rates(preNetNum), postnet_fwd_rates(postNetNum);
	preNet.GetNeuronType(prenet_types);
	postNet.GetNeuronType(postnet_types);
	for (int i = 0; i < preNetNum; i ++) {
		if (prenet_types[i]) {
			prenet_fwd_rates[i] = {prenet_rate_exc, 0.0};
		} else {
			prenet_fwd_rates[i] = {prenet_rate_inh, 0.0};
		}
	}
	for (int i = 0; i < postNetNum; i ++) {
		if (postnet_types[i]) {
			postnet_fwd_rates[i] = {postnet_rate_exc, 0.0};
		} else {
			postnet_fwd_rates[i] = {postnet_rate_inh, 0.0};
		}
	}
	if (preType) {
		preNet.InitializeExternalPoissonProcess(prenet_fwd_rates, maximum_time, atoi(m_map_config["PreNetExternalDrivingSeed"].c_str()));
	} else preNet.InitializeInternalPoissonRate(prenet_fwd_rates);
	if (postType) {
		postNet.InitializeExternalPoissonProcess(postnet_fwd_rates, maximum_time, atoi(m_map_config["PostNetExternalDrivingSeed"].c_str()));
	} else postNet.InitializeInternalPoissonRate(postnet_fwd_rates);
	// Set feedforward driving strength;
	preNet.SetF(true, atof(m_map_config["PreNetDrivingStrength"].c_str()));
	postNet.SetF(true, atof(m_map_config["PostNetDrivingStrength"].c_str()));

	// Set synaptic strength of presynaptic neural network;
	preNet.SetS(true, atof(m_map_config["PreNetSynapticStrengthExcitatory"].c_str()));
	postNet.SetS(true, atof(m_map_config["PostNetSynapticStrengthExcitatory"].c_str()));

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
	double sampling_rate = 1.0 / atof(m_map_config["RecordingTimingStep"].c_str());
	// Define the shape of data;
	size_t preNetshape[2], postNetshape[2];
	preNetshape[0] = tmax / dt;
	postNetshape[0] = tmax /dt;
	preNetshape[1] = preNetNum;
	postNetshape[1] = postNetNum;

	// Define file path for output data;
	FILEWRITE preV(dir + "preV.bin");
	FILEWRITE postV(dir + "postV.bin");
	FILEWRITE preI(dir + "preI.bin");
	FILEWRITE postI(dir + "postI.bin");
	// Initialize the shape of files:
	preV.SetSize(preNetshape);
	postV.SetSize(postNetshape);
	preI.SetSize(preNetshape);
	postI.SetSize(postNetshape);

	// char cr = (char)13;
	start = clock();
	int progress;
	while (t < tmax) {
		UpdateSystemState(preNet, postNet, conMat, t, dt);
		t += dt;
		// Output temporal data;
		if (floor(t * sampling_rate) == t * sampling_rate) {
			preNet.OutPotential(preV);
			preNet.OutCurrent(preI);
			postNet.OutPotential(postV);
			postNet.OutCurrent(postI);
		}
		if (floor(t * 10 / tmax) > progress) {
			progress = floor(t * 10/tmax);
			printf(">> Processing ... %d0%\n", progress);
		}
		// cout << cr;
	}

	// OUTPUTS:
	//preNet.SaveNeuron(dir + "preNeuron.csv");
	preNet.SaveConMat(dir + "preMat.csv");
	preNet.OutSpikeTrains(dir + "rasterPre.csv");
	//postNet.SaveNeuron(dir + "postNeuron.csv");
	postNet.SaveConMat(dir + "postMat.csv");
	postNet.OutSpikeTrains(dir + "rasterPost.csv");

	finish = clock();

	// COUNTS:
	printf(">> It takes %.2fs\n", (finish - start)*1.0 / CLOCKS_PER_SEC);
	return 0;
}
