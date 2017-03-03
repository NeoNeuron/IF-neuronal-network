//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-02-21 16:05:00
//	Description: test program for Class Neuron;
//***************
#include"neuron.h"
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
	double t = 0, dt = pow(2, -2), tmax = 200;
	double rateE = 0.5;
	double rateI = 0;
	double v;
	ofstream data;
	data.open("data.txt");
	for (int i = 0; i < 10; i++) {
		vector<double> x, y;
		GenerateExternalPoissonSequence(rateE, tmax, 1, x);
		GenerateExternalPoissonSequence(rateI, tmax, 2, y);
		cell.SetDrivingType(true);
		 // to be modified;
		while (t < tmax) {
			v = cell.dynamic(t, dt, x, y);
			t += dt;
			if (abs(floor(2 * t) - 2 * t) < 1e-15) {
				data << setprecision(20) << (double)v << "\t";
			}			
			//cout << '\b' << '\b' << '\b' << '\b' << '\b' << '\b' << '\b' << '\b' << '\b' << '\b';
			//cout << setprecision(4) << (double)(100 * t / tmax);			
		}
		//cell.outG();
		data << endl;
		cell.reset();
		dt /= 2;
		t = 0;
		cout << i << '\t';
	}
	cout << endl;
	finish = clock();
	// readout runing time;
	double ToTtime;
	ToTtime = (finish - start) / CLOCKS_PER_SEC;
	if (ToTtime < 60) {
		cout << "It takes " << (double)ToTtime << "s" << endl;
	} else if (ToTtime >= 60 && ToTtime < 3600) {
		int  MIN;
		double S;
		MIN = floor(ToTtime / 60);
		S = ToTtime - MIN * 60;
		cout << "It takes " << MIN << " min " << (double)S << " s" << endl;
	} else {
		int MIN, H;
		double S;
		H = ToTtime / 3600;
		MIN = (ToTtime - H * 3600) / 60;
		S = ToTtime - H * 3600 - MIN * 60;
		cout << "It takes " << H << " h " << MIN << " min " << (double)S << " s" << endl;
	}
	return 0;
}
