#include <cmath>
#include "../include/io.h"
#include "../include/spike.h"
using namespace std;

void Spike2Bool(vector<double>& spikes, vector<bool> & binary_spikes, double tmax, double dt) {
	int T = ceil(tmax / dt);
	binary_spikes.resize(T, false);
	for (vector<double>::iterator it = spikes.begin(); it != spikes.end(); it++) {
		int index = floor(*it / dt);
		if (index == T) index --;
		binary_spikes[index] = true;
	}
}

void Truncate(vector<double>& spikes, double* t_range) {
	vector<double>::iterator it = spikes.begin();
	while (it != spikes.end()) {
		if (*it <= t_range[0]) it = spikes.erase(it);
		else if (*it > t_range[1]) {
			spikes.erase(it, spikes.end());
			break;
		} else {
			*it -= t_range[0];
			it++;
		}
	}
}

void OutSpikeTrain(string path, vector<bool>& spikes) {
	Print1D(path, spikes, "trunc", 1);
}
