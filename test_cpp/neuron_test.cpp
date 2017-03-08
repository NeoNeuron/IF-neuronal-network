//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-02-21 16:05:00
//	Description: test program for Class Neuron;
//***************
#include"../neuron.h"
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
	double rateI = 0;
	double v;
	ofstream data;
	data.open("../data_new.txt");
	for (int i = 0; i < 15; i++) {
		vector<double> impt_e, impt_i;
		GenerateExternalPoissonSequence(rateE, tmax, 1, impt_e);
		GenerateExternalPoissonSequence(rateI, tmax, 2, impt_i);
		cell.SetDrivingType(true);
		while (t < tmax) {
			v = cell.UpdateNeuronalState(t, dt, impt_e, impt_i);
			t += dt;
			if (abs(floor(2 * t) - 2 * t) < 1e-15) {
				data << setprecision(20) << (double)v << "\t";
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
