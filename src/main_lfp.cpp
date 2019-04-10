// ***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2019-04-06
//	Description: main file for lfp.h and lfp.cpp
//***************
#include "../include/lfp.h"
#include "../include/io.h"
#include "../include/get-config.h"
#include "../include/common_header.h"
using namespace std;

//	Function of calculating LFP with point current source model in 1-D loop network case;
//
//	arguments:
//
//	argv[1] = dir of neural data file;
//	argv[2] = path of output LFP file;
//
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
	int neuron_num = atoi(m_map_config["NeuronNumber"].c_str());
	vector<int> list(neuron_num);
	for (int i = 0; i < neuron_num; i ++) list[i] = i;
	printf(">> %d connected neuron contribute to LFP\n", neuron_num);
	fflush(stdout);

	// Calculate the spatial weights
	vector<double> spatial_weights;
	CalculateSpatialWeight(m_map_config, spatial_weights);

	//	Choose objective time range;
	double t_range[2]; // t_range[0] = t_min; t_range[1] = t_max;
	t_range[0] = atof(m_map_config["TimeRangeMin"].c_str());
	t_range[1] = atof(m_map_config["TimeRangeMax"].c_str());
	printf(">> Time range = (%.2f, %.2f] ms\n", t_range[0], t_range[1]);
	fflush(stdout);

	// Processing LFP data;
	vector<double> lfp;
	string LFP_type = m_map_config["LFPType"];
	CalculateLFP(argv[1], lfp, list, LFP_type, spatial_weights, t_range, atof(m_map_config["SamplingTimingStep"].c_str()));

	//	Output lfp:
	Print1DBin(argv[2], lfp, "trunc");
	finish = clock();
	// counting time;
	printf("[-] LFP generation : %3.3f s\n", (finish - start)*1.0 / CLOCKS_PER_SEC);
	return 0;
}
