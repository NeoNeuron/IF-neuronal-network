#include "../include/io.h"
#include <iostream>
#include <sstream>
#include <fstream>
// #include <algorithm>
#include <numeric>
using namespace std;

double STA(vector<bool>& x, vector<double>& y) {
  double product = inner_product(x.begin(), x.end(), y.begin(), 0.0);
  return product / accumulate(x.begin(), x.end(), 0.0);
}

int main(int argc, const char* argv[]) {
  if (argc != 4) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
  // Loading data;
  string spike_path = argv[1];
  string lfp_path = argv[2];
  vector<bool> spike;
  vector<double> lfp;
  Read1DBin(spike_path, spike, 0, 0);
  Read1DBin(lfp_path, lfp, 0, 0);
  // Set time range;
  string buffer;
	istringstream range_in(argv[3]);
	getline(range_in, buffer, ',');
	int ntd = atoi(buffer.c_str());
	getline(range_in, buffer, ',');
	int ptd = atoi(buffer.c_str());
  // Calculate STA:
  vector<double> sta(ntd + ptd + 1);
  size_t res = ntd > ptd ? ntd : ptd;
	vector<bool> spike_copy(spike.begin(), spike.end() - res);
	vector<double> lfp_copy(lfp.begin(), lfp.end() - res);
  // No shift;
	sta[ntd] = STA(spike_copy, lfp_copy);
	// Negative shift;
	for (size_t i = 0; i < ntd; i++) {
		spike_copy.erase(spike_copy.begin());
		spike_copy.insert(spike_copy.end(), *(spike.end() - res + i));
		sta[ntd - i - 1] = STA(spike_copy, lfp_copy);
	}
	// Positive shift;
	spike_copy.clear();
	spike_copy.insert(spike_copy.end(), spike.begin(), spike.end() - res);
	for (size_t i = 0; i < ptd; i++) {
		lfp_copy.erase(lfp_copy.begin());
		lfp_copy.insert(lfp_copy.end(), *(lfp.end() - res + i));
		sta[ntd + i + 1] = STA(spike_copy, lfp_copy);
	}
  // Output STA:
  ofstream data_out;
	cout << ">> Outputing data ... " << endl;
	data_out.open("./data/sta/sta.csv");
	data_out << "timelag,sta" << endl;
	for (int i = 0; i < ntd + ptd + 1; i++) {
		data_out << i - ntd << ',' << setprecision(15) << (double)sta[i] << '\n';
	}
	data_out.close();

  finish = clock();
	// COUNTS:
	double ToTtime;
	ToTtime = (finish - start)*1.0 / CLOCKS_PER_SEC;
	printf(">> It takes %.2fs\n", ToTtime);
  return 0;
}
