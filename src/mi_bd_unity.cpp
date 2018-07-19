//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2018-01-27
//	Description: Mutual information analysis program; version 1.0
//***************
#include "../include/spike.h"
#include "../include/lfp.h"
#include "../include/mi_uniform.h"
#include "../include/io.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using namespace std;
int myrandom(int i) {return rand()%i;}
// compact function to calculate mutual information between multi-type signal
//	arguments:
//	argv[1] = path for bool series;
//	argv[2] = path for double series;
// 	argv[3] = index of spike;
// 	argv[4] = index of lfp;
// 	argv[5] = time range of series;
// 	argv[6] = bin size of series;
//	argv[7] = range of timelag;
//	argv[8] = bin size of pdf of the double variable;
int main(int argc, const char* argv[]) {
	if (argc != 9) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// Parameters:
	double trange[2];
	stringstream itrange(argv[5]);
	string buffer;
	for (size_t i = 0; i < 2; i++) {
		getline(itrange, buffer, ',');
		trange[i] = atof(buffer.c_str());
	}
	printf(">> Time range = (%.2f, %.2f] ms\n", trange[0], trange[1]);
	double dt = atof(argv[6]);
	// Set time range;
	size_t delay_range[2];
	istringstream i_delay_range(argv[7]);
	getline(i_delay_range, buffer, ',');
	int ntd = atoi(buffer.c_str());
	delay_range[0] = ntd;
	getline(i_delay_range, buffer, ',');
	int ptd = atoi(buffer.c_str());
	delay_range[1] = ptd;
	double binsize = atof(argv[8]);
	// INPUT NEURONAL DATA:
	vector<double> spikes;
	Read1D(argv[1], spikes, atoi(argv[3]), 0);
	// Truncate the spiking series;
	Truncate(spikes, trange);
	// Convert double to binary;
	vector<bool> binary_spikes;
	double tmax = trange[1] - trange[0];
	Spike2Bool(spikes, binary_spikes, tmax, dt);
	// LFP:
	vector<int> list;
	stringstream inlist(argv[4]);
	while (getline(inlist, buffer, ',')) {
		list.push_back(atoi(buffer.c_str()));
	}
	int neu_num = list.size();
	printf(">> %d connected neuron contribute to LFP\n", neu_num);
	cout << ">> Processing LFP ..." << endl;
	vector<double> lfp;
	LFP(argv[2], lfp, list, trange);

	// Calculate mutual information;
	cout << ">> Calculating TDMI ..." << endl;
	vector<double> tdmi;
	TDMI(binary_spikes, lfp, tdmi, delay_range, binsize);
	// shuffle flag;
	srand(unsigned(time(0)));
	random_shuffle(binary_spikes.begin(), binary_spikes.end(), myrandom);
	vector<double> tdmi_shuffle;
	TDMI(binary_spikes, lfp, tdmi_shuffle, delay_range, binsize);

	//	Output data:
	ofstream data_out;
	cout << ">> Outputing data ... " << endl;
	data_out.open("./data/mi/mi_bd_unity.csv");
	data_out << "timelag,mi,mi_shuffle" << endl;
	for (int i = 0; i < ntd + ptd + 1; i++) {
		data_out << i - ntd << ',' << setprecision(15) << (double)tdmi[i] << ',' << (double)tdmi_shuffle[i] << '\n';
	}
	data_out.close();

	finish = clock();
	// Time counting:
	cout << "[-] total TDMI calculation takes " << (finish - start)*1.0 / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
