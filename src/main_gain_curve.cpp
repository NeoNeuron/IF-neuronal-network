//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-09-02
//	Description: calculate the gain curve of the neuronal model;
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
	string net_config_path = "./doc/gain_curve.ini";
  map<string, string> m_map_config;
  ReadConfig(net_config_path,m_map_config);
  cout << ">> [Config.ini]:\n#####\n";
	PrintConfig(m_map_config);
	cout << "#####\n";
	double *dym_val;
	NeuronSim cell(m_map_config["NeuronType"], dym_val);
	double t = 0, dt = atof(m_map_config["TimingStep"].c_str()), tmax = atof(m_map_config["MaximumTime"].c_str());
	bool driving_type;
	istringstream(m_map_config["DrivingType"]) >> boolalpha >> driving_type;
	cell.SetDrivingType(driving_type);
	double rate_exc = atof(m_map_config["DrivingRate"].c_str());
	vector<double> in_E;
	if (driving_type) {
		GenerateExternalPoissonSequence(rate_exc, tmax, atoi(m_map_config["ExternalDrivingSeed"].c_str()), in_E);
	} else {
		srand(time(NULL));
		cell.SetPoissonRate(true, rate_exc);
		cell.SetPoissonRate(false, 0);
	}
	
	// Set driving strength;
	double s_min = atof(m_map_config["DrivingStrengthMin"].c_str());
	double s_max = atof(m_map_config["DrivingStrengthMax"].c_str());
	double ds = atof(m_map_config["ds"].c_str());
	vector<double> s_vec = {s_min};
	while (s_vec.back() + ds <= s_max) {
		s_vec.push_back(s_vec.back() + ds);
	} 

	start = clock();
	// SETUP DYNAMICS:
	//double recording_rate = 2;
	// Define the shape of data;
	//size_t shape[2];
	//shape[0] = tmax * recording_rate;
	//shape[1] = 1;
	// Define file path for output data;
	//string V_path = dir + "V.bin";
	//string I_path = dir + "I.bin";
	// Initialize files:
	//ofstream V_file, I_file;
	//V_file.open(V_path.c_str(), ios::binary);
	//I_file.open(I_path.c_str(), ios::binary);
	//V_file.write((char*)shape, 2*sizeof(size_t));
	//I_file.write((char*)shape, 2*sizeof(size_t));
	//double V, I;
	ofstream raster_file("./tmp/gain_curve.csv");
	raster_file.close();
	vector<double> new_spikes;
	vector<double> spike_train;
	vector<Spike> pe_spike;
	vector<Spike>::iterator it;
	for (int i = 0; i < s_vec.size(); i ++) {
		cell.Reset(dym_val);
		cell.SetDrivingType(true);
		pe_spike.resize(in_E.size());
		for (int j = 0; j < in_E.size(); j ++) {
			pe_spike[j].type = true;
			pe_spike[j].t = in_E[j];
			pe_spike[j].s = s_vec[i];
		}
		it = pe_spike.begin();
		while (t < tmax) {
			cell.UpdateNeuronalState(dym_val, t, dt, pe_spike, it, new_spikes);
			//if (new_spikes.size() > 0) cout << new_spikes.size() << '\t';
			if (new_spikes.size() > 0) cell.Fire(t, new_spikes);
			cell.CleanUsedInputs(dym_val, dym_val, t + dt);
			t += dt;
			//if (abs(recording_rate*t - floor(recording_rate*t)) == 0) {
			//cout << cell.GetPotential(dym_val) << '\t';
			//	//I = cell.OutTotalCurrent(dym_val);
			//	V_file.write((char*)&V, sizeof(double));
			//	//I_file.write((char*)&I, sizeof(double));
			//}
		}
		//cout << endl;
		//V_file.close();
		//I_file.close();

		//cell.GetCycle();
		//cout<< endl;
		// OUTPUTS:
		cell.OutSpikeTrain(spike_train);
		Print1D("./tmp/gain_curve.csv", spike_train, "app", 0);
		printf(">> s = %.2e, mean firing rate: %.2f Hz\n", s_vec[i], (double)spike_train.size()*1000.0/tmax);
		t = 0;
	}
	finish = clock();
	// readout runing time;
	cout << "It takes " << (finish - start)*1.0 / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
