#include <cmath>
#include "../include/io.h"
#include "../include/spike.h"
using namespace std;

void Spike2Bool(vector<double>& spikes, vector<bool> & binary_spikes, double tmax, double dt) {
	int T = ceil(tmax / dt);
	binary_spikes.resize(T, false);
	for (vector<double>::iterator it = spikes.begin(); it != spikes.end(); it ++) {
		int index = floor(*it / dt);
		if (index == T) index --;
		binary_spikes[index] = true;
	}
}

void Truncate(vector<double>& spikes, double* range) {
	vector<double> spikes_temp;
	for (vector<double>::iterator it = spikes.begin(); it < spikes.end(); it ++) {
		if (*it > range[0] && *it <= range[1]) spikes_temp.push_back(*it - range[0]);
		else if (*it > range[1]) break;
	}
	spikes.clear();
	spikes = spikes_temp;
}
