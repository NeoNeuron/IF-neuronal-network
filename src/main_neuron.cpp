//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-09-02
//	Description: simulating program for Class Neuron;
//***************
#include "../include/neuron.h"
#include "../include/io.h"
#include "../include/get-config.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>
#include <stdexcept>
using namespace std;

int main(int argc, const char* argv[]) {
	if (argc != 2) throw runtime_error("wrong number of args");
	clock_t start, finish;
	// 	Setup directory for output files;
	//	it must be existing dir;
	string dir;
	dir = argv[1];

	// Loading config.ini:
	string net_config_path = "./doc/config_neuron.ini";
  map<string, string> m_map_config;
  ReadConfig(net_config_path,m_map_config);
  cout << ">> [Config.ini]:\n#####\n";
	PrintConfig(m_map_config);
	cout << "#####\n";
	Neuron cell;
	double t = 0, dt = atof(m_map_config["TimingStep"].c_str()),  tmax = atof(m_map_config["MaximumTime"].c_str());
	bool neuron_type, driving_type;
	istringstream(m_map_config["Type"]) >> boolalpha >> neuron_type;
	istringstream(m_map_config["DrivingType"]) >> boolalpha >> driving_type;
	cell.SetNeuronType(neuron_type);
	cell.SetDrivingType(driving_type);
	double rate_exc = atof(m_map_config["DrivingRate"].c_str());
	vector<double> in_E, in_I;
	if (driving_type) {
		GenerateExternalPoissonSequence(rate_exc, tmax, 1, in_E);
		GenerateExternalPoissonSequence(0, tmax, 2, in_I);
	} else {
		srand(time(NULL));
		cell.SetPoissonRate(true, rate_exc);
	}
	// Set driving strength;
	cell.SetFeedforwardStrength(true, atof(m_map_config["DrivingStrength"].c_str()));

	start = clock();
	// SETUP DYNAMICS:
	double recording_rate = 2;
	// Define the shape of data;
	size_t shape[2];
	shape[0] = tmax * recording_rate;
	shape[1] = 1;
	// Define file path for output data;
	string V_path = dir + "V.bin";
	string I_path = dir + "I.bin";
	// Initialize files:
	ofstream V_file, I_file;
	V_file.open(V_path.c_str(), ios::binary);
	I_file.open(I_path.c_str(), ios::binary);
	V_file.write((char*)shape, 2*sizeof(size_t));
	I_file.write((char*)shape, 2*sizeof(size_t));
	double V, I;
	while (t < tmax) {
		cell.UpdateNeuronalState(t, dt, in_E, in_I);
		t += dt;
		if (abs(recording_rate*t - floor(recording_rate*t)) == 0) {
			V = cell.GetPotential();
			I = cell.OutTotalCurrent();
			V_file.write((char*)&V, sizeof(double));
			I_file.write((char*)&I, sizeof(double));
		}
	}
	V_file.close();
	I_file.close();

	// OUTPUTS:
	vector<double> spike_train;
	string raster_path = dir + "raster.csv";
	cell.OutSpikeTrain(spike_train);
	Print1D(raster_path, spike_train, "trunc", 0);

	finish = clock();
	// readout runing time;
	cout << "It takes " << (finish - start)*1.0 / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
