//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-06-03
//	Description: Mutual information analysis program; version 1.0
//***************
#include "../include/mi.h"
#include "../include/io.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using namespace std;
// compact function to calculate mutual information between multi-type signal
//	arguments:
//	argv[1] = mode code;
//		0: spike to spike
//		1: spike to lfp
//		2: lfp to lfp
//	argv[2] = timing step for TDMI;
//	argv[3] = delay time range;
int main(int argc, const char* argv[]) {
	if (argc != 4) {
		throw runtime_error("wrong number of args");
	}
	clock_t start, finish;
	start = clock();
	// Preparing input args;
	double dt = atof(argv[2]);
	string range = argv[3];
	string::size_type pos = range.find_first_of(',', 0 );
	int negative_time_delay = atoi(range.substr(0, pos).c_str());
	range.erase(0, pos + 1);
	int positive_time_delay = atoi(range.c_str());
	range = "";
	printf(">> dt = %f ms\n", dt);
	printf(">> maximum negative time delay = %d\n", negative_time_delay);
	printf(">> maximum positive time delay = %d\n", positive_time_delay);
	// Judge the running mode:
	if (argv[1] == "0") {
		vector<double> raster_x, raster_y;
		string path;
		path = "./data/raster/raster_x.csv";
		Read1D(path, 0, 1, raster_x);
		path = "./data/raster/raster_y.csv";
		Read1D(path, 0, 1, raster_y);

		double tmax;
		if (raster_x.back() > raster_y.back()) tmax = ceil(raster_x.back());
		else tmax = ceil(raster_y.back());
		cout << ">> Calculating ordered TDMI ... " << endl;
		vector<double> tdmi;
		TDMI(raster_x, raster_y, dt, tmax, negative_time_delay, positive_time_delay, tdmi);
		//	Output data:
		ofstream data_out;
		cout << ">> Outputing data ... " << endl;
		data_out.open("./data/mi/mi_ss.csv");
		data_out << "timelag,mi" << endl;
		for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i++) {
			data_out << (double)dt*(i - negative_time_delay) << ',' << (double)tdmi[1];
			if (i < positive_time_delay + negative_time_delay) data_out << endl;
		}
		data_out.close();
	} else if (argv[1] == "1") {
		// INPUT NEURONAL DATA:
		vector<double> raster, lfp;

		// DATA OF PRELAYER NEURON:
		string path;
		path = "./data/raster/raster.csv";
		Read1D(path, 0, 1, raster);
		path = "./data/lfp/lfp.csv";
		Read1D(path, 0, 1, lfp);

		double sampling_dt = 0.03125;
		cout << ">> Calculating ordered TDMI ... " << endl;
		vector<double> tdmi_ordered;
		TDMI(raster, lfp, dt, sampling_dt, negative_time_delay, positive_time_delay, tdmi_ordered, false);
		cout << ">> Calculating swapped TDMI ... " << endl;
		vector<double> tdmi_random;
		TDMI(raster, lfp, dt, sampling_dt, negative_time_delay, positive_time_delay, tdmi_random, true);

		//	Output data:
		ofstream data_out;
		cout << ">> Outputing data ... " << endl;
		data_out.open("./data/mi/mi_sl.csv");
		data_out << "timelag,ordered,random" << endl;
		for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i++) {
			data_out << (double)dt*(i - negative_time_delay) << ',' << (double)tdmi_ordered[i] << ',' << (double)tdmi_random[i];
			if (i < positive_time_delay + negative_time_delay) data_out << endl;
		}
		data_out.close();
	} else if (argv[1] == "2") {
		// INPUT NEURONAL DATA:
		vector<double> lfp_x, lfp_y;

		// DATA OF PRELAYER NEURON:
		string path;
		path = "./data/lfp/lfp_x.csv";
		Read1D(path, 0, 1, lfp_x);
		path = "./data/lfp/lfp_y.csv";
		Read1D(path, 0, 1, lfp_y);

		double sampling_dt = 0.03125;
		// take average;
		int np = floor(dt/sampling_dt);
		int nd = floor(lfp_x.size()/np);
		vector<double> lfp_x_ave(nd, 0), lfp_y_ave(nd, 0);
		for (int i = 0; i < nd; i++) {
			for (int j = 0; j < np; j++) {
				lfp_x_ave[i] += lfp_x[i*np + j]/np;
				lfp_y_ave[i] += lfp_y[i*np + j]/np;
			}
		}
		cout << ">> Calculating ordered TDMI ... " << endl;
		vector<double> tdmi;
		TDMI_adaptive(lfp_x_ave, lfp_y_ave, negative_time_delay, positive_time_delay, tdmi);
		//	Output data:
		ofstream data_out;
		cout << ">> Outputing data ... " << endl;
		data_out.open("./data/mi/mi_ll.csv");
		data_out << "timelag,mi" << endl;
		for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i++) {
			data_out << (double)dt*(i - negative_time_delay) << ',' << (double)tdmi[i];
			if (i < positive_time_delay + negative_time_delay) data_out << endl;
		}
		data_out.close();
	}

	finish = clock();
	// Time counting:
	double ToTtime;
	ToTtime = (finish - start) / CLOCKS_PER_SEC;
	cout << "It takes " << (double)ToTtime << "s" << endl;
	return 0;
}
