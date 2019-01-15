//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-09-02
//	Description: simulating program for Class Neuron;
//***************
#include "../include/neuron.h"
#include "../include/io.h"
#include "../include/get-config.h"
#include "../include/poisson_generator.h"
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

	// ========================================
	// Loading config.ini:
	// ========================================
	string net_config_path = "./doc/config_neuron.ini";
  map<string, string> m_map_config;
  ReadConfig(net_config_path,m_map_config);
  cout << ">> [Config.ini]:\n#####\n";
	PrintConfig(m_map_config);
	cout << "#####\n";
	// ========================================
	// Prepare config parameters:
	// ========================================
	string neuron_type = m_map_config["NeuronType"];
	double t = 0;
	double dt = atof(m_map_config["TimingStep"].c_str());
	double tmax = atof(m_map_config["MaximumTime"].c_str());
	int poisson_num = atoi(m_map_config["PoissonNum"].c_str());
	double const_I = atof(m_map_config["ConstDrive"].c_str());

	// ========================================
	// Initialize Neuron and Poisson generator;
	// ========================================
	double *dym_val, *dym_val_new;
	dym_val_new = new double[4];
	NeuronSim cell(neuron_type, dym_val);
	cell.SetConstDrive(const_I);
	vector<PoissonGenerator> pgs(poisson_num);
	for (int i = 0; i < 4; i ++) dym_val_new[i] = dym_val[i];
	double pr, ps;
	for (int i = 0; i < poisson_num; i ++) {
		pr = atof( m_map_config["pr" + to_string(i)].c_str() );
		ps = atof( m_map_config["ps" + to_string(i)].c_str() );
		pgs[i].SetRate(pr);
		pgs[i].SetStrength(ps);
		pgs[i].SetOuput( m_map_config["PoissonDir"] + "pg" + to_string(i) + ".csv" );
	}
	// ========================================
	// Generate Poisson:
	// ========================================
	srand( atoi(m_map_config["pSeed"].c_str()) );
	//srand((int)time(NULL));
	vector<Spike> pe_spike;
	for (int i = 0; i < poisson_num; i ++) {
		pgs[i].GenerateNewPoisson( tmax, pe_spike );
	}
	// erase the spike at t = 0;
	pe_spike.erase(pe_spike.begin(), pe_spike.begin() + poisson_num);
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
	//I_file.open(I_path.c_str(), ios::binary);
	V_file.write((char*)shape, 2*sizeof(size_t));
	//I_file.write((char*)shape, 2*sizeof(size_t));
	double V, I;
	vector<double> new_spikes;
	vector<Spike>::iterator it = pe_spike.begin();

	// ========================================
	// Start simulation; 
	// ========================================
	start = clock();
	while (t < tmax) {
		cell.UpdateNeuronalState(dym_val_new, t, dt, pe_spike, it, new_spikes);
		if (new_spikes.size() > 0) cell.Fire(t, new_spikes);
		cell.CleanUsedInputs(dym_val, dym_val_new, t + dt);
		t += dt;
		if (abs(recording_rate*t - floor(recording_rate*t)) == 0) {
			V = cell.GetPotential(dym_val);
			//cout << V << ',';
			//I = cell.OutTotalCurrent(dym_val);
			V_file.write((char*)&V, sizeof(double));
			//I_file.write((char*)&I, sizeof(double));
		}
	}
	V_file.close();
	//I_file.close();

	cell.GetCycle();
	cout<< endl;
	// OUTPUTS:
	vector<double> spike_train;
	cell.OutSpikeTrain(spike_train);
	cout << "Mean firing rate: " << spike_train.size() * 1000.0 / tmax << endl;
	Print1D(dir + "raster.csv", spike_train, "trunc", 0);

	finish = clock();
	// readout runing time;
	cout << "It takes " << (finish - start)*1.0 / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
