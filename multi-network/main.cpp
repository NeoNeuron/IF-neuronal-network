//*************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-08 17:21:36
//	Description: test program for multi-network simulation;
//*************************

#include "multi_network.h"
#include "get-config.h"
#include <cmath>
#include <ctime>
#include <sstream>

using namespace std;

int main() {
	clock_t start, finish;

	// 	Setup directory for output files;
	//	it must be existing dir;
	string dir;
	cout << "Enter existing directory of output files: " << endl;
	cin >> dir;
	
	start = clock();
	// SETUP NETWORKS:

	string m_sPath = "./config.ini";  
  map<string, string> m_mapConfig;  
  ReadConfig(m_sPath,m_mapConfig);  
  PrintConfig(m_mapConfig);
	
	int preNetNum, postNetNum, preNetDensity, postNetDensity;
	preNetNum = atoi(m_mapConfig["PreNetNeuronNumber"].c_str());
	postNetNum = atoi(m_mapConfig["PostNetNeuronNumber"].c_str());
	preNetDensity = atoi(m_mapConfig["PreNetConnectingDensity"].c_str());
	postNetDensity = atoi(m_mapConfig["PostNetConnectingDensity"].c_str());
	NeuronalNetwork preNet(preNetNum, preNetDensity), postNet(postNetNum, postNetDensity);
	double maximum_time = atof(m_mapConfig["MaximumTime"].c_str());
	//cout << "Input the maximum simulating time (ms):" << endl;
	//cin >> maximum_time;
	string neuronFileName, conMatFileName;
	//neuronFileName = dir + "preNeuron.txt";
	//conMatFileName = dir + "preMat.txt";
	//preNet.load(neuronFileName, conMatFileName);
	//neuronFileName = dir + "postNeuron.txt";
	//conMatFileName = dir + "postMat.txt";
	//postNet.load(neuronFileName, conMatFileName);
	preNet.InitializeNeuronalType(atof(m_mapConfig["PreNetTypeProbability"].c_str()), atoi(m_mapConfig["PreNetTypeSeed"].c_str()));
	postNet.InitializeNeuronalType(atof(m_mapConfig["PostNetTypeProbability"].c_str()), atoi(m_mapConfig["PostNetTypeSeed"].c_str()));

	bool type;
	istringstream(m_mapConfig["PreNetDrivingType"]) >> boolalpha >> type;
	preNet.SetDrivingType(type);
	preNet.InitializeExternalPoissonProcess(true, atof(m_mapConfig["PreNetDrivingRate"].c_str()), maximum_time, atoi(m_mapConfig["PreNetExternalDrivingSeed"].c_str()));
	preNet.InitializeExternalPoissonProcess(false, 0, maximum_time, 0);

	istringstream(m_mapConfig["PostNetDrivingType"]) >> boolalpha >> type;
	postNet.SetDrivingType(type);
	postNet.InitializeExternalPoissonProcess(true, atof(m_mapConfig["PostNetDrivingRate"].c_str()), maximum_time, atoi(m_mapConfig["PostNetExternalDrivingSeed"].c_str()));
	postNet.InitializeExternalPoissonProcess(false, 0, maximum_time, 0);

	// SETUP CONNECTIVITY MAT
	vector<vector<bool> > conMat;
	double conP = atof(m_mapConfig["ConnectingProbability"].c_str());
	srand(atoi(m_mapConfig["ConnectingSeed"].c_str()));
	vector<bool> addRow;
	ofstream conMatFile;
	string conMatFilename = dir + "conMat.txt";
	conMatFile.open(conMatFilename.c_str());
	for (int i = 0; i < preNetNum; i++) {
		addRow.clear();
		for (int j = 0; j < postNetNum; j++) {
			if (rand() / (RAND_MAX*1.0) < conP) {
				addRow.push_back(true);
				conMatFile << true << "\t";
			} else {
				addRow.push_back(false);
				conMatFile << false << "\t";
			}
		}
		conMat.push_back(addRow);
		conMatFile << endl;
	}
	conMatFile.close();

	// SETUP DYNAMICS:
	double t = 0, dt = atof(m_mapConfig["TimingStep"].c_str()), tmax = maximum_time;
	ofstream preV, preGE, preGI, postV, postGE, postGI;
	string preV_filename = dir + "preV.txt";
	string preGE_filename = dir + "preGE.txt";
	string preGI_filename = dir + "preGI.txt";
	string postV_filename = dir + "postV.txt";
	string postGE_filename = dir + "postGE.txt";
	string postGI_filename = dir + "postGI.txt";
	preV.open(preV_filename.c_str());
	preGE.open(preGE_filename.c_str());
	preGI.open(preGI_filename.c_str());
	postV.open(postV_filename.c_str());	
	postGE.open(postGE_filename.c_str());
	postGI.open(postGI_filename.c_str());

	vector<vector<double> > readme;
	char cr = (char)13;
	double progress;
	cout << "Begin ... " << endl;
	while (t < tmax) {
		UpdateSystemState(preNet, postNet, conMat, t, dt);
		t += dt;
		readme.clear();
		preNet.OutputTemporalParameters(readme);
		for (int i = 0; i < preNetNum; i++) {
			preV << setprecision(15) << (double)readme[0][i] << '\t';
			preGE << setprecision(15) << (double)readme[1][i] << '\t';
			preGI << setprecision(15) << (double)readme[2][i] << '\t';
		}
		preV << endl;
		preGE << endl;
		preGI << endl;

		readme.clear();
		postNet.OutputTemporalParameters(readme);
		for (int i = 0; i < postNetNum; i++) {
			postV << setprecision(15) << (double)readme[0][i] << '\t';
			postGE << setprecision(15) << (double)readme[1][i] << '\t';
			postGI << setprecision(15) << (double)readme[2][i] << '\t';
		}
		postV << endl;
		postGE << endl;
		postGI << endl;
		progress = t * 100.0 / tmax;
		cout << cr << ">> Processing ...";
		printf("%6.2f", progress);
		cout << "%";
	}
	cout << endl << "END!" << endl;

	preV.close();
	preGE.close();
	postGE.close();
	postV.close();
	preGI.close();
	postGI.close();

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
	preNet.OutputSpikeTrains(spikes);
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
	postNet.OutputSpikeTrains(spikes);
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
		cout << "It takes " << (double)ToTtime << "s" << endl;
	} else if (ToTtime >= 60 && ToTtime < 3600) {
		int  MIN;
		double S;
		MIN = floor(ToTtime / 60);
		S = ToTtime - MIN * 60;
		cout << "It takes " << MIN << " min " << (double)S << " s" << endl;
	} else {
		int MIN, H;
		double S;
		H = ToTtime / 3600;
		MIN = (ToTtime - H * 3600) / 60;
		S = ToTtime - H * 3600 - MIN * 60;
		cout << "It takes " << H << " h " << MIN << " min " << (double)S << " s" << endl;
	}
	return 0;
}
