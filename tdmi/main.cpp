//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-08 22:36:03
//	Description: Mutual information analysis program; version 1.0
//***************
#include "mi.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cmath>

using namespace std;

int main() {
	clock_t start, finish;
	

	// INPUT NEURONAL DATA:
	ifstream data;
	vector<double> raster_test;
	vector<double> LFP_test;
	string file_dir = "../lfp/file-txt/";

	// DATA OF PRELAYER NEURON:

	string file_name;
	file_name = file_dir + "lfp_test.txt";
	ReadData(file_name, LFP_test);
	file_name = file_dir + "raster_test.txt";
	ReadData(file_name, raster_test);

	// TIME DELAY

	
	int expected_occupancy, negative_time_delay, positive_time_delay;
	double dt, dtSampling;
	cout << ">> Input expected occupancy: ";
	cin >> expected_occupancy;
	cout << ">> Input dt: ";
	cin >> dt;
	cout << ">> Input maximum negative time delay: ";
	cin >> negative_time_delay;
	cout << ">> Input maximum positive time delay: ";
	cin >> positive_time_delay;

	double sampling_dt = 0.03125;

	start = clock();
	cout << ">> Calculate TDMI based on ordered spike train ... " << endl;	
	vector<double> tdmi_ordered;
	TDMI(raster_test, LFP_test, expected_occupancy, dt, sampling_dt, negative_time_delay, positive_time_delay, tdmi_ordered, false);
	cout << ">> Calculate TDMI based on swapped spike train ... " << endl;	
	vector<double> tdmi_random;
	TDMI(raster_test, LFP_test, expected_occupancy, dt, sampling_dt, negative_time_delay, positive_time_delay, tdmi_random, true);

	//	Output data:
	ofstream data_out;
	cout << ">> Outputing data ... " << endl;
	data_out.open("./file-dat/tdmi_ordered.dat");
	for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i++) {
		data_out << (double)dt*(i - negative_time_delay) << "\t" << (double)tdmi_ordered[i];
		if (i < positive_time_delay + negative_time_delay) data_out << endl;
	}
	data_out.close();
	data_out.open("./file-dat/tdmi_rand.dat");
	for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i++) {
		data_out << (double)dt*(i - negative_time_delay) << "\t" << (double)tdmi_random[i];
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