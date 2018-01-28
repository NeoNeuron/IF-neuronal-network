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

//	Simulation program for two network system;
//	arguments:
//	argv[1] = Outputing directory for neural data;
int main(int argc, const char* argv[]) {
	if (argc != 2) {
		throw runtime_error("wrong number of args");
	}
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
	vector<vector<double> > prenet_fwd_rates(preNetNum), postnet_fwd_rates(postNetNum);
	preNet.GetNeuronType(prenet_types);
	postNet.GetNeuronType(postnet_types);
	for (int i = 0; i < preNetNum; i ++) {
		if (prenet_types[i]) {
			prenet_fwd_rates[i].push_back(prenet_rate_exc);
			prenet_fwd_rates[i].push_back(0.0);
		} else {
			prenet_fwd_rates[i].push_back(prenet_rate_inh);
			prenet_fwd_rates[i].push_back(0.0);
		}
	}
	for (int i = 0; i < postNetNum; i ++) {
		if (postnet_types[i]) {
			postnet_fwd_rates[i].push_back(postnet_rate_exc);
			postnet_fwd_rates[i].push_back(0.0);
		} else {
			postnet_fwd_rates[i].push_back(postnet_rate_inh);
			postnet_fwd_rates[i].push_back(0.0);
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
	// Define file path for output data;
	// string preV_path = dir + "preV.csv";
	// string preGE_path = dir + "preGE.csv";
	// string preGI_path = dir + "preGI.csv";
	// string postV_path = dir + "postV.csv";
	// string postGE_path = dir + "postGE.csv";
	// string postGI_path = dir + "postGI.csv";
	string preI_path = dir + "preI.csv";
	string postI_path = dir + "postI.csv";
	// string preEI_path = dir + "preEI.csv";
	// string preII_path = dir + "preII.csv";
	// string postEI_path = dir + "postEI.csv";
	// string postII_path = dir + "postII.csv";
	// Initialize files:
	// ofstream preV, postV;
	ofstream preI, postI;
	// preV.open(preV_path.c_str());
	// preV.close();
	// preGE.open(preGE_path.c_str());
	// preGE.close();
	// preGI.open(preGI_path.c_str());
	// postGE.close();
	// postV.open(postV_path.c_str());
	// postV.close();
	// postGE.open(postGE_path.c_str());
	// postGE.close();
	// postGI.open(postGI_path.c_str());
	// postGI.close();
	preI.open(preI_path.c_str());
	preI.close();
	postI.open(postI_path.c_str());
	postI.close();
	// preEI.open(preEI_path.c_str());
	// preEI.close();
	// preII.open(preII_path.c_str());
	// preII.close();
	// postEI.open(postEI_path.c_str());
	// postEI.close();
	// postII.open(postII_path.c_str());
	// postII.close();

	// char cr = (char)13;
	// double progress;
	while (t < tmax) {
		UpdateSystemState(preNet, postNet, conMat, t, dt);
		t += dt;
		// Output temporal data;
		// preNet.OutPotential(preV_path);
		preNet.OutCurrent(preI_path);
		// preNet.OutPartialCurrent(preEI, true);
		// preNet.OutPartialCurrent(preII, false);
		// postNet.OutPotential(postV_path);
		postNet.OutCurrent(postI_path);
		// postNet.OutPartialCurrent(postEI, true);
		// postNet.OutPartialCurrent(postII, false);

		// progress = t * 100.0 / tmax;
		// cout << cr;
		// printf(">> Processing ... %6.2f", progress);
		// cout << "%";
	}
	cout << endl;
	//
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
