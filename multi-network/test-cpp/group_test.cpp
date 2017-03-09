//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-02-21 16:06:04
//	Description: test program for Class NeuronalNetwork;
//***************
#include"../group.h"
#include<cstdio>
#include<ctime>
#include<cmath>

using namespace std;

int main() {
	clock_t start, finish;
	start = clock();
	NeuronalNetwork cells(20, 2);
	double t = 0, dt = pow(2, -1), tmax = 2000;
	double rateE = 1.5;
	double rateI = 0;
	cells.InitializeNeuronalType(0.5, 1);
	vector<double> x;
	ofstream data;
	for (int i = 0; i < 10; i++) {
		cells.RestoreNeurons();
		string fileName;
		char ss[3];
		sprintf(ss, "%d", i + 3);
		fileName += "data_new";
		fileName += ss;
		fileName += ".txt";
		cells.InitializeExternalPoissonProcess(true, rateE, tmax, 1);
		cells.InitializeExternalPoissonProcess(false, rateI, tmax, 4);
		cells.SetDrivingType(true); // external type;
		//cells.InitializeInternalPoissonRate(false, rateI);
		const char* char_filename = fileName.c_str();
		data.open(char_filename);
		while (t < tmax) {
			cells.UpdateNetworkState(t, dt);
			t += dt;
			x.clear();
			cells.OutputPotential(x);
			if (abs(floor(2 * t) - 2 * t) < 1e-15) {
				for (vector<double>::iterator it = x.begin(); it != x.end(); it++) {
					data << setprecision(20) << (double)*it << '\t';
				}
				data << endl;
			}
		}
		data.close();
		cout << endl;
		cout << i << endl;
		t = 0;
		dt /= 2;
	}	
	vector<vector<double> > spikes;
	cells.OutputSpikeTrains(spikes);
	data.open("raster.txt");
	for (vector<vector<double> >::iterator it = spikes.begin(); it != spikes.end(); it++) {
		for (vector<double>::iterator itt = it->begin(); itt != it->end(); itt++) {
			data << setprecision(5) << (double)*itt << '\t';
		}
		data << endl;
	}
	data.close();

	vector<bool> types;
	cells.OutputNeuronType(types);
	data.open("type.txt");
	for (vector<bool>::iterator it = types.begin(); it != types.end(); it++) {
		data << *it << endl;
	}
	data.close();

	finish = clock();
	cout << endl;
	double ToTtime;
	ToTtime = (finish - start)*1.0 / CLOCKS_PER_SEC;
	cout << "The program was excecuted within " << (double)ToTtime << " s. " << endl;	
	return 0;
}
