// ***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2019-04-06
//	Description: main file for lfp.h and lfp.cpp
//***************
#include "../include/lfp.h"
#include "../include/io.h"
#include "../include/common_header.h"
using namespace std;

//	Function of calculating LFP with point current source model in 1-D loop network case;
//
//	arguments:
//
//	argv[1] = path of neural data file;
//	argv[2] = path of output LFP file;
//	argv[3] = index(indices) of target neurons (seperated by comma);
//	argv[4] = time range;
//	argv[5] = sampling time step;
//
int main(int argc, const char* argv[]) {
	if (argc != 6) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// Analyze listing series;
	vector<int> list;
	stringstream inlist(argv[3]);
	string buffer;
	while (getline(inlist, buffer, ',')) {
		list.push_back(atoi(buffer.c_str()));
	}
	int neuron_num = list.size();
	printf(">> %d connected neuron contribute to LFP\n", neuron_num);
	fflush(stdout);

	// Define uniform spatial weights
	vector<double> spatial_weights(neuron_num, 1.0);

	//	Choose objective time range;
	double t_range[2]; // t_range[0] = t_min; t_range[1] = t_max;
	stringstream strange(argv[4]);
	for (int i = 0; i < 2; i ++) {
		getline(strange, buffer, ',');
		t_range[i] = atof(buffer.c_str());
	}
	printf(">> Time range = (%.2f, %.2f] ms\n", t_range[0], t_range[1]);
	fflush(stdout);

	// Processing LFP data;
	vector<double> lfp;
	string LFP_type = "tot";
	LFP(argv[1], lfp, list, LFP_type, spatial_weights, t_range, atof(argv[5]));

	//	Output lfp:
	Print1DBin(argv[2], lfp, "trunc");
	finish = clock();
	// counting time;
	printf("[-] LFP generation : %3.3f s\n", (finish - start)*1.0 / CLOCKS_PER_SEC);
	fflush(stdout);

	return 0;
}
