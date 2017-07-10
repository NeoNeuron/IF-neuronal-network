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

using namespace std;
//	Function of calculating time-delayed mutual information with point;
//	arguments:
//	argv[1] = timing step for TDMI;
//	argv[2] = lower bond of time-delay range;
//	argv[3] = upper bond of time-delay range;
int main(int argc, const char* argv[]) {
	if (argc != 4) {
		throw runtime_error("wrong number of args");
	}
	clock_t start, finish;
	start = clock();
	// INPUT NEURONAL DATA:
	ifstream data;
	vector<double> raster;
	vector<double> lfp;
	string file_dir = "./lfp/file-txt/";

	// DATA OF PRELAYER NEURON:
	string path;
	path = file_dir + "raster.txt";
	Read1D(path, 0, 1, raster);
	path = file_dir + "lfp.txt";
	Read1D(path, 0, 1, lfp);
	// Preparing input args;
	double dt = atof(argv[1]);
	int negative_time_delay = atoi(argv[2]);
	int positive_time_delay = atoi(argv[3]);

	printf(">> dt = %f ms\n", dt);
	printf(">> maximum negative time delay = %d\n", negative_time_delay);
	printf(">> maximum positive time delay = %d\n", positive_time_delay);


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
	data_out.open("./tdmi/file-txt/tdmi_ordered.txt");
	for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i++) {
		data_out << (double)dt*(i - negative_time_delay) << ',' << (double)tdmi_ordered[i];
		if (i < positive_time_delay + negative_time_delay) data_out << endl;
	}
	data_out.close();

	data_out.open("./tdmi/file-txt/tdmi_rand.txt");
	for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i++) {
		data_out << (double)dt*(i - negative_time_delay) << ',' << (double)tdmi_random[i];
		if (i < positive_time_delay + negative_time_delay) data_out << endl;
	}
	data_out.close();

	finish = clock();
	// Time counting:
	double ToTtime;
	ToTtime = (finish - start) / CLOCKS_PER_SEC;
	cout << "It takes " << (double)ToTtime << "s" << endl;
	return 0;
}
