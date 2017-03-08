//*************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-08 17:21:36
//	Description: test program for multi-network simulation;
//*************************

#include "multi_network.h"
#include <cmath>
#include <ctime>

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
	
	int preNetNum, postNetNum;
	preNetNum = 100;
	postNetNum = 100;
	NeuronalNetwork preNet(preNetNum, 3), postNet(postNetNum, 3);
	double maximum_time = 10000;	
	//cout << "Input the maximum simulating time (ms):" << endl;
	//cin >> maximum_time;
	string neuronFileName, conMatFileName;
	//neuronFileName = dir + "preNeuron.txt";
	//conMatFileName = dir + "preMat.txt";
	//preNet.load(neuronFileName, conMatFileName);
	//neuronFileName = dir + "postNeuron.txt";
	//conMatFileName = dir + "postMat.txt";
	//postNet.load(neuronFileName, conMatFileName);
	preNet.InitializeNeuronalType(0.8, 1);
	postNet.InitializeNeuronalType(0.8, 2);

	preNet.SetDrivingType(true);
	preNet.InitializeExternalPoissonProcess(true, 1.5, maximum_time, 5);
	preNet.InitializeExternalPoissonProcess(false, 0, maximum_time, 5);

	postNet.SetDrivingType(true);
	postNet.InitializeExternalPoissonProcess(true, 1.5, maximum_time, 6);
	postNet.InitializeExternalPoissonProcess(false, 0, maximum_time, 6);

	// SETUP CONNECTIVITY MAT
	vector<vector<bool> > conMat;
	double conP = 0.1;
	srand(100);
	vector<bool> addRow;
	ofstream conMatFile;
	string conMatFilename = dir + "conMat.txt";
	const char* char_connectivity_matrix_filename = conMatFilename.c_str();
	conMatFile.open(char_connectivity_matrix_filename);
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
	double t = 0, dt = pow(2, -5), tmax = maximum_time;
	ofstream preV, preGE, preGI, postV, postGE, postGI;
	string preV_filename = dir + "preV.txt";
	string preGE_filename = dir + "preGE.txt";
	string preGI_filename = dir + "preGI.txt";
	string postV_filename = dir + "postV.txt";
	string postGE_filename = dir + "postGE.txt";
	string postGI_filename = dir + "postGI.txt";
	const char* char_preV_filename = preV_filename.c_str();
	const char* char_preGE_filename = preGE_filename.c_str();
	const char* char_preGI_filename = preGI_filename.c_str();
	const char* char_postV_filename = postV_filename.c_str();
	const char* char_postGE_filename = postGE_filename.c_str();
	const char* char_postGI_filename = postGI_filename.c_str();	
	preV.open(char_preV_filename);
	preGE.open(char_preGE_filename);
	postGE.open(char_postGE_filename);
	postV.open(char_postV_filename);
	preGI.open(char_preGI_filename);
	postGI.open(char_postGI_filename);
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
	const char* char_pre_raster_filename = pre_raster_filename.c_str();
	const char* char_post_raster_filename = post_raster_filename.c_str();
	spikes.clear();
	preNet.OutputSpikeTrains(spikes);
	data.open(char_pre_raster_filename);
	for (vector<vector<double> >::iterator it = spikes.begin(); it != spikes.end(); it++) {
		for (vector<double>::iterator itt = it->begin(); itt != it->end(); itt++) {
			data << setprecision(5) << (double)*itt << '\t';
		}
		data << endl;
	}
	data.close();

	spikes.clear();
	postNet.OutputSpikeTrains(spikes);
	data.open(char_post_raster_filename);
	for (vector<vector<double> >::iterator it = spikes.begin(); it != spikes.end(); it++) {
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
