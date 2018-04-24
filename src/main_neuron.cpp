//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-09-02
//	Description: simulating program for Class Neuron;
//***************
#include "../include/neuron.h"
#include "../include/io.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <stdexcept>
using namespace std;

int main(int argc, const char* argv[]) {
	if (argc != 3) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	Neuron cell;
	double t = 0, dt = pow(2, -5), tmax = atof(argv[1]);
	srand(time(NULL));
	int ntrials = atoi(argv[2]);
	double rateI = 0;
	vector<double> I;
	ofstream ofile;
	ofile.open("./data/currents.csv", ios::trunc);
	ofile.close();
	for (int i = 0; i < ntrials; i++) {
		cell.SetDrivingType(false);
		cell.SetPoissonRate(1.5);
		I.clear();
		while (t < tmax) {
			cell.UpdateNeuronalState(t, dt);
			I.push_back(cell.OutTotalCurrent());
			t += dt;
		}
		cell.Reset();
		t = 0;
		cout << i << ' ';
		Print1D("./data/currents.csv", I, "app", 0);
	}
	cout << endl;
	finish = clock();
	// readout runing time;
	double ToTtime;
	cout << "It takes " << (finish - start)*1.0 / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
