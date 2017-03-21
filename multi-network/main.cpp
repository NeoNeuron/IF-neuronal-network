//*************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-13 15:07:52
//	Description: test program for multi-network simulation;
//*************************

#include "multi_network.h"
#include "get-config.h"
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <sstream>

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
	string net_config_path = "./multi-network/config.ini";
  map<string, string> m_map_config;  
  ReadConfig(net_config_path,m_map_config);
  cout << ">> [Config.ini]:" << endl;
	PrintConfig(m_map_config);  
	cout << endl;
	
	int preNetNum, postNetNum, preNetDensity, postNetDensity;
	preNetNum = atoi(m_map_config["PreNetNeuronNumber"].c_str());
	postNetNum = atoi(m_map_config["PostNetNeuronNumber"].c_str());
	preNetDensity = atoi(m_map_config["PreNetConnectingDensity"].c_str());
	postNetDensity = atoi(m_map_config["PostNetConnectingDensity"].c_str());
	NeuronalNetwork preNet(preNetNum, preNetDensity), postNet(postNetNum, postNetDensity);
	double maximum_time = atof(m_map_config["MaximumTime"].c_str());
	//cout << "Input the maximum simulating time (ms):" << endl;
	//cin >> maximum_time;
	string neuronFileName, conMatFileName;
	//neuronFileName = dir + "preNeuron.txt";
	//conMatFileName = dir + "preMat.txt";
	//preNet.load(neuronFileName, conMatFileName);
	//neuronFileName = dir + "postNeuron.txt";
	//conMatFileName = dir + "postMat.txt";
	//postNet.load(neuronFileName, conMatFileName);
	preNet.InitializeNeuronalType(atof(m_map_config["PreNetTypeProbability"].c_str()), atoi(m_map_config["PreNetTypeSeed"].c_str()));
	cout << "in pre-network." << endl;
	postNet.InitializeNeuronalType(atof(m_map_config["PostNetTypeProbability"].c_str()), atoi(m_map_config["PostNetTypeSeed"].c_str()));
	cout << "in post-network." << endl;

	bool type;
	istringstream(m_map_config["PreNetDrivingType"]) >> boolalpha >> type;
	preNet.SetDrivingType(type);
	preNet.InitializeExternalPoissonProcess(true, atof(m_map_config["PreNetDrivingRate"].c_str()), maximum_time, atoi(m_map_config["PreNetExternalDrivingSeed"].c_str()));
	preNet.InitializeExternalPoissonProcess(false, 0, maximum_time, 0);

	istringstream(m_map_config["PostNetDrivingType"]) >> boolalpha >> type;
	postNet.SetDrivingType(type);
	postNet.InitializeExternalPoissonProcess(true, atof(m_map_config["PostNetDrivingRate"].c_str()), maximum_time, atoi(m_map_config["PostNetExternalDrivingSeed"].c_str()));
	postNet.InitializeExternalPoissonProcess(false, 0, maximum_time, 0);

	// SETUP CONNECTIVITY MAT
	vector<vector<bool> > conMat(preNetNum);
	double conP = atof(m_map_config["ConnectingProbability"].c_str());
	srand(atoi(m_map_config["ConnectingSeed"].c_str()));
	ofstream conMatFile;
	string conMatFilename = dir + "conMat.txt";
	conMatFile.open(conMatFilename.c_str());
	for (int i = 0; i < preNetNum; i++) {
		conMat[i].resize(postNetNum);
		for (int j = 0; j < postNetNum; j++) {
			if (rand() / (RAND_MAX*1.0) < conP) {
				conMat[i][j] = true;
				conMatFile << true << "\t";
			} else {
				conMat[i][j] = false;
				conMatFile << false << "\t";
			}
		}
		conMatFile << endl;
	}
	conMatFile.close();

	// SETUP DYNAMICS:
	double t = 0, dt = atof(m_map_config["TimingStep"].c_str()), tmax = maximum_time;
	ofstream preV, preGE, preGI, postV, postGE, postGI;
	ofstream preI, postI;
	string preV_filename = dir + "preV.txt";
	string preGE_filename = dir + "preGE.txt";
	string preGI_filename = dir + "preGI.txt";
	string postV_filename = dir + "postV.txt";
	string postGE_filename = dir + "postGE.txt";
	string postGI_filename = dir + "postGI.txt";
	string preI_filename = dir + "preI.txt";
	string postI_filename = dir + "postI.txt";
	preV.open(preV_filename.c_str());
	preGE.open(preGE_filename.c_str());
	preGI.open(preGI_filename.c_str());
	postV.open(postV_filename.c_str());	
	postGE.open(postGE_filename.c_str());
	postGI.open(postGI_filename.c_str());
	preI.open(preI_filename.c_str());
	postI.open(postI_filename.c_str());

	vector<vector<double> > readme;
	vector<double> current;
	char cr = (char)13;
	double progress;
	while (t < tmax) {
		UpdateSystemState(preNet, postNet, conMat, t, dt);
		t += dt;
		readme.clear();
		current.clear();
		preNet.OutTemporalParameters(readme);
		preNet.OutCurrent(current);
		for (int i = 0; i < preNetNum; i++) {
			preV << setprecision(15) << (double)readme[i][0] << '\t';
			preGE << setprecision(15) << (double)readme[i][1] << '\t';
			preGI << setprecision(15) << (double)readme[i][2] << '\t';
			preI << setprecision(15) << (double)current[i] << '\t';
		}
		preV << endl;
		preGE << endl;
		preGI << endl;
		preI << endl;

		readme.clear();
		current.clear();
		postNet.OutTemporalParameters(readme);
		postNet.OutCurrent(current);
		for (int i = 0; i < postNetNum; i++) {
			postV << setprecision(15) << (double)readme[i][0] << '\t';
			postGE << setprecision(15) << (double)readme[i][1] << '\t';
			postGI << setprecision(15) << (double)readme[i][2] << '\t';
			postI << setprecision(15) << (double)current[i] << '\t';
		}
		postV << endl;
		postGE << endl;
		postGI << endl;
		postI << endl;
		progress = t * 100.0 / tmax;
		cout << cr;
		printf(">> Processing ... %6.2f", progress);
		cout << "%";
	}
	cout << endl;

	preV.close();
	preGE.close();
	postGE.close();
	postV.close();
	preGI.close();
	postGI.close();
	preI.close();
	postI.close();

	neuronFileName = dir + "preNeuron.txt";
	conMatFileName = dir + "preMat.txt";
	preNet.Save(neuronFileName, conMatFileName);
	neuronFileName = dir + "postNeuron.txt";
	conMatFileName = dir + "postMat.txt";
	postNet.Save(neuronFileName, conMatFileName);

	// OUTPUTS:
	vector<vector<double> > spikes;
	ofstream data;
	string pre_raster_filename = dir + "rasterPre.txt";
	string post_raster_filename = dir + "rasterPost.txt";
	spikes.clear();
	preNet.OutSpikeTrains(spikes);
	data.open(pre_raster_filename.c_str());
	int neuron_indices = 0;
	for (vector<vector<double> >::iterator it = spikes.begin(); it != spikes.end(); it++) {
		data << (int) neuron_indices << '\t';
		neuron_indices ++;
		for (vector<double>::iterator itt = it->begin(); itt != it->end(); itt++) {
			data << setprecision(5) << (double)*itt << '\t';
		}
		data << endl;
	}
	data.close();

	spikes.clear();
	postNet.OutSpikeTrains(spikes);
	data.open(post_raster_filename.c_str());
	neuron_indices = 0;
	for (vector<vector<double> >::iterator it = spikes.begin(); it != spikes.end(); it++) {
		data << (int) neuron_indices << '\t';
		neuron_indices ++;
		for (vector<double>::iterator itt = it->begin(); itt != it->end(); itt++) {
			data << setprecision(5) << (double)*itt << '\t';
		}
		data << endl;
	}
	data.close();

	finish = clock();

	// COUNTS:
	double ToTtime;
	ToTtime = (finish - start) / CLOCKS_PER_SEC;
	if (ToTtime < 60) {
		printf(">> It takes %.2fs\n", ToTtime);
	} else if (ToTtime >= 60 && ToTtime < 3600) {
		int  MIN;
		double S;
		MIN = floor(ToTtime / 60);
		S = ToTtime - MIN * 60;
		printf(">> It takes %dmin %.2fs\n", MIN, S);
	} else {
		int MIN, H;
		double S;
		H = ToTtime / 3600;
		MIN = (ToTtime - H * 3600) / 60;
		S = ToTtime - H * 3600 - MIN * 60;
		printf(">> It takes %dh %dmin %.2fs\n", H, MIN, S);
	}
	cout << ">> END!" << endl;
	return 0;
}
