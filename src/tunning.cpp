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
	if (argc != 5) {
		throw runtime_error("wrong number of args");
	}
	// 	Setup directory for output files;
	//	it must be existing dir;
	double df = atof(argv[1]);
	double duf = atof(argv[2]);
	string frange = argv[3];
	string::size_type pos = frange.find_first_of(',', 0 );
	double f_min = atof(frange.substr(0, pos).c_str());
	frange.erase(0, pos + 1);
	double f_max = atof(frange.c_str());
	frange = "";
	string ufrange = argv[4];
	pos = ufrange.find_first_of(',', 0 );
	double uf_min = atof(ufrange.substr(0, pos).c_str());
	ufrange.erase(0, pos + 1);
	double uf_max = atof(ufrange.c_str());
	ufrange = "";

	Neuron cell;
	double t = 0, dt = 0.1, tmax = 10000;
	double f = f_min, uf = uf_min;
	vector<double> impt_e, impt_i;
	vector<double> spike_train;
	vector<vector<double> > firing_rates;
	vector<double> newline;
	double u;
	while (f <= f_max) {
		newline.clear();
		while (uf <= uf_max) {
			u = uf / f;
			// cout << f << endl;
			cell.SetDrivingType(false);
			cell.SetPoissonRate(true, u);
			cell.SetFeedforwardConductance(true, f);
			// cout << u << ',' << f << endl;
			while (t < tmax) {
				cell.UpdateNeuronalState(t, dt, impt_e, impt_i);
				// cout << v << endl;
				t += dt;
			}
			cell.OutSpikeTrain(spike_train);
			// cout << spike_train.size()*1000.0/tmax << endl;
			newline.push_back(spike_train.size()*1000.0/tmax);
			cell.Reset();
			t = 0;
			uf += duf;
		}
		uf = uf_min;
		f += df;
		firing_rates.push_back(newline);
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
	Print2D(path, firing_rates, "trunc");

	return 0;
}
