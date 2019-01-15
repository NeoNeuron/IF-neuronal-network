//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-02-21 16:05:00
//	Description: test program for Class Neuron;
//***************
#include "../include/neuron.h"
#include "../include/io.h"
#include "../include/poisson_generator.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

int main() {
	clock_t start, finish;
	start = clock();
	double *dym_val; 
	double *dym_val_new;
	dym_val_new = new double[4];
	NeuronSim cell("LIF_G", dym_val);
	double t = 0, dt = 0.5, tmax = 4000;
	double rateE = 1.5;
	double v;
	ofstream data;
	data.open("./tmp/data_new.txt");
	PoissonGenerator pg;
	pg.SetRate(rateE);
	pg.SetStrength(9e-3);
	vector<Spike> pe_spike;
	pg.GenerateNewPoisson(tmax, pe_spike);
	vector<Spike>::iterator it_cur;
	double sampling_rate = 2;
	vector<double> new_spikes;
	vector<double> spike_train;
	for (int i = 0; i < 8; i++) {
		cell.Reset(dym_val);
		for (int j = 0; j < 4; j ++) dym_val_new[j] = dym_val[j];
		it_cur = pe_spike.begin();
		while (t < tmax) {
			cell.UpdateNeuronalState(dym_val_new, t, dt, pe_spike, it_cur, new_spikes);
			if (new_spikes.size() > 0) cell.Fire(t, new_spikes);
			v = cell.CleanUsedInputs(dym_val, dym_val_new, t + dt);
			t += dt;
			if (abs(floor(sampling_rate * t) - sampling_rate * t) < 1e-15) {
				data << setprecision(20) << (double)v << ",";
			}
		}
		data << endl;
		cell.OutSpikeTrain(spike_train);
		Print1D("./tmp/data_neuron_raster.csv", spike_train, "app", 0);
		printf(">> dt = %.2e, mean firing rate: %.2f Hz\n",dt,(double)spike_train.size()*1000.0/tmax);
		dt /= 2;
		t = 0;
	}
	finish = clock();
	// readout runing time;
	cout << "It takes " << (finish - start) / CLOCKS_PER_SEC << " s" << endl;
	return 0;
}
