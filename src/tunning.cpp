//*************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-13 15:07:52
//	Description: test program for multi-network simulation;
//*************************

#include "../include/neuron.h"
// #include "../include/get-config.h"
#include "../include/io.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <sstream>
#include <stdexcept>
using namespace std;

//	Simulation program for single network system;
//	arguments:
//	argv[1] = Outputing directory for neural data;
int main(int argc, const char* argv[]) {
	if (argc != 3) {
		throw runtime_error("wrong number of args");
	}
	// 	Setup directory for output files;
	//	it must be existing dir;
	double df = atof(argv[1]);
	string frange = argv[2];
	string::size_type pos = frange.find_first_of(',', 0 );
	double f_min = atof(frange.substr(0, pos).c_str());
	frange.erase(0, pos + 1);
	double f_max = atof(frange.c_str());
	frange = "";

	Neuron cell;
	double t = 0, dt = 0.1, tmax = 10000;
	double f = f_min;
	vector<double> impt_e, impt_i;
	vector<double> spike_train;
	vector<double> firing_rates;
	double v;
	while (f <= f_max) {
		// cout << f << endl;
		cell.SetDrivingType(false);
		cell.SetPoissonRate(true, 7.5);
		cell.SetPoissonRate(false, 0);
		cell.SetFeedforwardConductance(true, f);
		while (t < tmax) {
			v = cell.UpdateNeuronalState(t, dt, impt_e, impt_i);
			// cout << v << endl;
			t += dt;
		}
		cell.OutSpikeTrain(spike_train);
		// cout << spike_train.size()*1000.0/tmax << endl;
		firing_rates.push_back(spike_train.size()*1000.0/tmax);
		cell.Reset();
		t = 0;
		f += df;
	}

	// char cr = (char)13;
	// double progress;
	// while (t < tmax) {
	// 	net.UpdateNetworkState(t, dt);
	// 	t += dt;
	// 	// Output temporal data;
	// 	net.OutPotential(V_path);
	// 	net.OutCurrent(I_path);
	//
	// 	progress = t * 100.0 / tmax;
	// 	cout << cr;
	// 	printf(">> Processing ... %6.2f", progress);
	// 	cout << "%";
	// }
	// cout << endl;

	// OUTPUTS:
	string path = "./data/tuninng.csv";
	Print1D(path, "trunc", 1, firing_rates);

	return 0;
}
