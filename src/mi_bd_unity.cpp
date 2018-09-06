//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2018-01-27
//	Description: Compact Mutual information analysis program; version 1.0
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
//
// arguments:
//
// argv[1] = directory of data files;
// argv[2] = index of spike;
// argv[3] = index of lfp;
// argv[4] = time range of series;
// argv[5] = bin size of series;
// argv[6] = range of timelag;
// argv[7] = bin size of pdf of the double variable;

int main(int argc, const char* argv[]) {
	if (argc != 8) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// Prepare data path;
	string dir = argv[1];
	string spike_file = dir + "raster.csv";
	string lfp_file = dir + "I.bin";
	string mi_file = dir + "mi_bd.csv";
	// Parameters:
	double trange[2];
	stringstream itrange(argv[4]);
	string buffer;
	for (size_t i = 0; i < 2; i++) {
		getline(itrange, buffer, ',');
		trange[i] = atof(buffer.c_str());
	}
	printf(">> Time range = (%.2f, %.2f] ms\n", trange[0], trange[1]);
	double dt = atof(argv[5]);
	// Set time range;
	size_t delay_range[2];
	istringstream i_delay_range(argv[6]);
	getline(i_delay_range, buffer, ',');
	int ntd = atoi(buffer.c_str());
	delay_range[0] = ntd;
	getline(i_delay_range, buffer, ',');
	int ptd = atoi(buffer.c_str());
	delay_range[1] = ptd;
	
	// Input neuronal data:
	
	clock_t clock_begin, clock_end;
	// Spike:
	clock_begin = clock(); 
	vector<double> spikes;
	Read1D(spike_file, spikes, atoi(argv[2]), 0);
	// Truncate the spiking series;
	Truncate(spikes, trange);
	clock_end = clock();
	// print time;
	cout << "[-] spike processing took " << (clock_end - clock_begin)*1.0 / CLOCKS_PER_SEC << "s" << endl;

	// Convert double to binary;
	vector<bool> binary_spikes;
	double tmax = trange[1] - trange[0];
	Spike2Bool(spikes, binary_spikes, tmax, dt);

	// LFP:
	clock_begin = clock(); 
	vector<int> list;
	stringstream inlist(argv[3]);
	while (getline(inlist, buffer, ',')) {
		list.push_back(atoi(buffer.c_str()));
	}
	int neu_num = list.size();
	printf(">> %d connected neuron contribute to LFP ...\n", neu_num);
	vector<double> lfp;
	LFP(lfp_file, lfp, list, trange);
	clock_end = clock();
	// print time;
	cout << "[-] LFP processing took " << (clock_end - clock_begin)*1.0 / CLOCKS_PER_SEC << "s" << endl;

	// Calculate mutual information;
	clock_begin = clock(); 
	// for mi_bd.out, param is binsize; mi_bd_2bins.out param is threshold;
	double param = atof(argv[7]);
	vector<double> tdmi;
  //TDMI(binary_spikes, lfp, tdmi, delay_range, param);
	TDMI2bins(binary_spikes, lfp, tdmi, delay_range, param);
	// shuffle flag;
	srand(unsigned(time(0)));
	random_shuffle(binary_spikes.begin(), binary_spikes.end(), myrandom);
	vector<double> tdmi_shuffle;
	//TDMI(binary_spikes, lfp, tdmi_shuffle, delay_range, param);
	TDMI2bins(binary_spikes, lfp, tdmi_shuffle, delay_range, param);
	clock_end = clock();
	// print time;
	cout << "[-] LFP processing took " << (clock_end - clock_begin)*1.0 / CLOCKS_PER_SEC << "s" << endl;

	//	Output data:
	ofstream data_out;
	data_out.open(mi_file.c_str());
	data_out << "#timelag,mi,mi_shuffle" << endl;
	for (int i = 0; i < ntd + ptd + 1; i++) {
		data_out << i - ntd << ',' << setprecision(15) << (double)tdmi[i] << ',' << (double)tdmi_shuffle[i] << '\n';
	}
	data_out.close();

	finish = clock();
	// Time counting:
	cout << "[-] total TDMI calculation takes " << (finish - start)*1.0 / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
