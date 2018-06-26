//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-02-21 16:05:00
//	Description: test program for Class Neuron;
//***************
#include"../include/neuron.h"
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include<ctime>
#include<cmath>

using namespace std;

int main() {
	clock_t start, finish;
	start = clock();
	Neuron cell;
	double t = 0, dt = 0.5, tmax = 2000;
	double rateE = 2;
	double rateI = 0.0;
	double v;
	ofstream data;
	data.open("./tmp/data_new.txt");
	vector<double> in_E, in_I;
	GenerateExternalPoissonSequence(rateE, tmax, 1, in_E);
	//GenerateExternalPoissonSequence(rateI, tmax, 2, in_I);
	vector<double> in_E_tmp, in_I_tmp;
	double sampling_rate = 2;
	double t_spike;
	for (int i = 0; i < 8; i++) {
		cell.Reset();
		in_E_tmp = in_E;
		in_I_tmp = in_I;
		cell.SetDrivingType(true);
		while (t < tmax) {
			//v = cell.UpdateNeuronalState(t, dt, in_E_tmp, in_I_tmp);
			t_spike = cell.TemporallyUpdateNeuronalState(t, dt, in_E_tmp, in_I_tmp);
			if (t_spike >= 0) {
				//cout << t_spike << ',';
				cell.Fire(t, t_spike);
				v = cell.UpdateNeuronalState(t + t_spike, dt - t_spike, in_E_tmp, in_I_tmp);
			} else {
				v = cell.UpdateNeuronalState(t, dt, in_E_tmp, in_I_tmp);
			}
			t += dt;
			if (abs(floor(sampling_rate * t) - sampling_rate * t) < 1e-15) {
				data << setprecision(20) << (double)v << ",";
			}
		}
		data << endl;
		dt /= 2;
		t = 0;
		cout << i << '\n';
		cout << endl;
	}
	cout << endl;
	vector<double> spikes;
	cell.OutSpikeTrain(spikes);
	cout << "Mean firing rate: " << (double)spikes.size()*1000.0/tmax << endl;
	finish = clock();
	// readout runing time;
	cout << "It takes " << (finish - start) / CLOCKS_PER_SEC << " s" << endl;
	return 0;
}
