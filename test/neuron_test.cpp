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
	double t = 0, dt = pow(2, -1), tmax = 2000;
	double rateE = 1.5;
	double v;
	ofstream data;
	data.open("./tmp/data_new.txt");
	vector<double> in;
	GenerateExternalPoissonSequence(rateE, tmax, 1, in);
	vector<double> in_tmp;
	int sampling_rate = 2;
	for (int i = 0; i < 8; i++) {
		in_tmp = in;
		cell.SetDrivingType(true);
		while (t < tmax) {
			v = cell.UpdateNeuronalState(t, dt, in_tmp);
			t += dt;
			if (abs(floor(sampling_rate * t) - sampling_rate * t) < 1e-15) {
				data << setprecision(20) << (double)v << ",";
			}
			//cout << setprecision(4) << (double)(100 * t / tmax);
		}
		data << endl;
		cell.Reset();
		dt /= 2;
		t = 0;
		cout << i << '\t';
	}
	cout << endl;
	finish = clock();
	// readout runing time;
	double ToTtime;
	ToTtime = (finish - start) / CLOCKS_PER_SEC;
	cout << "It takes " << (double)ToTtime << "s" << endl;
	return 0;
}
