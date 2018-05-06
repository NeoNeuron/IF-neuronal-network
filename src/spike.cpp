#include <cmath>
#include "../include/io.h"
#include "../include/spike.h"
using namespace std;

void Spike2Bool(vector<double>& spikes, vector<bool> & binary_spikes, double tmax, double dt) {
	size_t T = ceil(tmax / dt);
	binary_spikes.resize(T, false);
	size_t index;
	for (vector<double>::iterator it = spikes.begin(); it != spikes.end(); it ++) {
		index = floor((*it + 0.5*dt) / dt) - 1;
		if (index >= 0) {
			if (index == T) index --;
			binary_spikes[index] = true;
		}
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
