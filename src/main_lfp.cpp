// ***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2018-09-02
//	Description: main file for lfp.h and lfp.cpp
//***************
#include "../include/lfp.h"
#include "../include/io.h"
#include "../include/get-config.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
using namespace std;

//	Function of calculating LFP with point current source model in 1-D loop network case;
//
//	arguments:
//
//	argv[1] = path of neural data file;
//	argv[2] = path of output LFP file;
//
double L2(vector<double>& a, vector<double>& b) {
	return sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1]));
}

int main(int argc, const char* argv[]) {
	if (argc != 3) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// Load config file:
	string net_config_path = "./doc/config_lfp.ini";
  map<string, string> m_map_config;
  ReadConfig(net_config_path,m_map_config);
  cout << ">> [Config.ini]:\n#####\n";
	PrintConfig(m_map_config);
	cout << "#####\n";
	
	// Analyze listing series;
	vector<int> list;
	int select_mode = atoi(m_map_config["SelectMode"].c_str());
	if (select_mode == 0) {
		stringstream inlist(m_map_config["SelectedList"].c_str());
		string buffer;
		while (getline(inlist, buffer, ',')) {
			list.push_back(atoi(buffer.c_str()));
		}
	} else if (select_mode == 1) {
		double select_p = atof(m_map_config["SelectProbability"].c_str());
		int neu_num = atoi(m_map_config["NeuronNumber"].c_str());
		if (select_p == 1) {
			list.resize(neu_num);
			for (int i = 0; i < neu_num; i ++) list[i] = i;
		} else {
			srand(atof(m_map_config["SelectSeed"].c_str()));
			double x;
			for (int i = 0; i < neu_num; i ++) {
				x = rand()/(RAND_MAX * 1.0);
				if (x < select_p) list.push_back(i);
			}
		}
	}
	int neuron_num = list.size();
	printf(">> %d connected neuron contribute to LFP\n", neuron_num);

	// Calculate the spatial weights of selected neurons:
	vector<double> electrode_pos = {atof(m_map_config["PosX"].c_str()), atof(m_map_config["PosY"].c_str())};
	int decay_order = atoi(m_map_config["DecayOrder"].c_str());
	// read the coordinate file;
	vector<vector<double> > coordinates;
	Read2D(m_map_config["CoorPath"], coordinates);
	vector<double> spatial_weights(neuron_num);
	double distance;
	for (int i = 0; i < neuron_num; i ++) {
		distance = L2(electrode_pos, coordinates[list[i]]);
		if (decay_order == 1) spatial_weights[i] = 1.0 / distance;
		else if (decay_order == 2) spatial_weights[i] = 1.0 / distance / distance;
		else {
			cout << "Not proper decay order\n";
			break;
		}
	}

	//spatial_weights.clear();
	//spatial_weights.resize(neuron_num, 1.0);

	//	Choose objective time range;
	double t_range[2]; // t_range[0] = t_min; t_range[1] = t_max;
	t_range[0] = atof(m_map_config["TimeRangeMin"].c_str());
	t_range[1] = atof(m_map_config["TimeRangeMax"].c_str());
	printf(">> Time range = (%.2f, %.2f] ms\n", t_range[0], t_range[1]);

	// Processing LFP data;
	vector<double> lfp;
	LFP(argv[1], lfp, list, spatial_weights, t_range, atof(m_map_config["SamplingTimingStep"].c_str()));

	//	Output lfp:
	Print1DBin(argv[2], lfp, "trunc");
	finish = clock();
	// counting time;
	cout << "[-] LFP generation takes " << (finish - start)*1.0 / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
