//*************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-13 15:07:52
//	Description: test program for multi-network simulation;
//*************************

#include "../include/group.h"
#include "../include/get-config.h"
#include "../include/io.h"
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
  cout << ">> [Config.ini]:" << endl;
	PrintConfig(m_map_config);
	cout << endl;

	int preNetNum, postNetNum, preNetDensity, postNetDensity;
	double pre_rewiring_probability, post_rewiring_probability;
	int pre_rewiring_seed, post_rewiring_seed;
	preNetNum = atoi(m_map_config["PreNetNeuronNumber"].c_str());
	postNetNum = atoi(m_map_config["PostNetNeuronNumber"].c_str());
	preNetDensity = atoi(m_map_config["PreNetConnectingDensity"].c_str());
	postNetDensity = atoi(m_map_config["PostNetConnectingDensity"].c_str());
	pre_rewiring_probability = atof(m_map_config["PreNetRewiringProbability"].c_str());
	post_rewiring_probability = atof(m_map_config["PostNetRewiringProbability"].c_str());
	pre_rewiring_seed = atoi(m_map_config["PreNetRewiringSeed"].c_str());
	post_rewiring_seed = atoi(m_map_config["PostNetRewiringSeed"].c_str());
	// Generate networks;
	NeuronalNetwork preNet(preNetNum), postNet(postNetNum);
	preNet.SetConnectingDensity(preNetDensity);
	postNet.SetConnectingDensity(postNetDensity);
	preNet.Rewire(pre_rewiring_probability, pre_rewiring_seed, true);
	postNet.Rewire(post_rewiring_probability, post_rewiring_seed, true);
	double maximum_time = atof(m_map_config["MaximumTime"].c_str());

	preNet.InitializeNeuronalType(atof(m_map_config["PreNetTypeProbability"].c_str()), atoi(m_map_config["PreNetTypeSeed"].c_str()));
	cout << "in pre-network." << endl;
	postNet.InitializeNeuronalType(atof(m_map_config["PostNetTypeProbability"].c_str()), atoi(m_map_config["PostNetTypeSeed"].c_str()));
	cout << "in post-network." << endl;

	bool type;
	istringstream(m_map_config["PreNetDrivingType"]) >> boolalpha >> type;
	preNet.SetDrivingType(type);
	preNet.InitializeExternalPoissonProcess(true, atof(m_map_config["PreNetDrivingRateExcitatory"].c_str()), atof(m_map_config["PreNetDrivingRateInhibitory"].c_str()), maximum_time, atoi(m_map_config["PreNetExternalDrivingSeed"].c_str()));
	preNet.InitializeExternalPoissonProcess(false, 0, 0, maximum_time, 0);
	preNet.SetF(true, atof(m_map_config["PreNetDrivingStrength"].c_str()));

	istringstream(m_map_config["PostNetDrivingType"]) >> boolalpha >> type;
	postNet.SetDrivingType(type);
	postNet.InitializeExternalPoissonProcess(true, atof(m_map_config["PostNetDrivingRateExcitatory"].c_str()), atof(m_map_config["PostNetDrivingRateInhibitory"].c_str()), maximum_time, atoi(m_map_config["PostNetExternalDrivingSeed"].c_str()));
	postNet.InitializeExternalPoissonProcess(false, 0, 0, maximum_time, 0);
	postNet.SetF(true, atof(m_map_config["PostNetDrivingStrength"].c_str()));

	// SETUP CONNECTIVITY MAT
	vector<vector<bool> > conMat(preNetNum);
	double conP = atof(m_map_config["ConnectingProbability"].c_str());
	srand(atoi(m_map_config["ConnectingSeed"].c_str()));
	for (int i = 0; i < preNetNum; i++) {
		conMat[i].resize(postNetNum);
		for (int j = 0; j < postNetNum; j++) {
			if (rand() / (RAND_MAX*1.0) < conP) conMat[i][j] = true;
			else conMat[i][j] = false;
		}
	}
	string conMat_path = dir + "conMat.txt";
	Print2D(conMat_path, "trunc", conMat);

	// SETUP DYNAMICS:
	double t = 0, dt = atof(m_map_config["TimingStep"].c_str()), tmax = maximum_time;
	// Define file path for output data;
	string preV_path = dir + "preV.txt";
	// string preGE_path = dir + "preGE.txt";
	// string preGI_path = dir + "preGI.txt";
	string postV_path = dir + "postV.txt";
	// string postGE_path = dir + "postGE.txt";
	// string postGI_path = dir + "postGI.txt";
	string preI_path = dir + "preI.txt";
	string postI_path = dir + "postI.txt";
	// string preEI_path = dir + "preEI.txt";
	// string preII_path = dir + "preII.txt";
	// string postEI_path = dir + "postEI.txt";
	// string postII_path = dir + "postII.txt";
	// Initialize files:
	ofstream preV, postV, preI, postI;
	preV.open(preV_path.c_str());
	preV.close();
	// preGE.open(preGE_path.c_str());
	// preGE.close();
	// preGI.open(preGI_path.c_str());
	// postGE.close();
	postV.open(postV_path.c_str());
	postV.close();
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

	char cr = (char)13;
	double progress;
	while (t < tmax) {
		UpdateSystemState(preNet, postNet, conMat, t, dt);
		t += dt;
		// Output temporal data;
		preNet.OutPotential(preV_path);
		preNet.OutCurrent(preI_path);
		// preNet.OutPartialCurrent(preEI, true);
		// preNet.OutPartialCurrent(preII, false);
		postNet.OutPotential(postV_path);
		postNet.OutCurrent(postI_path);
		// postNet.OutPartialCurrent(postEI, true);
		// postNet.OutPartialCurrent(postII, false);

		progress = t * 100.0 / tmax;
		cout << cr;
		printf(">> Processing ... %6.2f", progress);
		cout << "%";
	}
	cout << endl;

	string preNeuron_path, preMat_path;
	preNeuron_path = dir + "preNeuron.txt";
	preMat_path = dir + "preMat.txt";
	preNet.Save(preNeuron_path, preMat_path);
	string postNeuron_path, postMat_path;
	postNeuron_path = dir + "postNeuron.txt";
	postMat_path = dir + "postMat.txt";
	postNet.Save(postNeuron_path, postMat_path);

	// OUTPUTS:
	string pre_raster_path = dir + "rasterPre.txt";
	string post_raster_path = dir + "rasterPost.txt";
	preNet.OutSpikeTrains(pre_raster_path);
	postNet.OutSpikeTrains(post_raster_path);

	finish = clock();

	// COUNTS:
	double ToTtime;
	ToTtime = (finish - start) / CLOCKS_PER_SEC;
	printf(">> It takes %.2fs\n", ToTtime);
	return 0;
}
