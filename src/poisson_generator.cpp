#include "../include/poisson_generator.h"
#include <cmath>
#include <cstdlib>

void PoissonGenerator::GenerateNewPoisson( double tmax, vector<Spike>& synaptic_driven) {
	Spike new_spike;
	new_spike.type = true;
	new_spike.s = strength_;
	double x, tLast = last_poisson_time_;
	while (tLast < tmax) {
		new_spike.t = tLast;
		synaptic_driven.push_back(new_spike);
		if (output_flag_) {
			outfile_ << setprecision(18) << tLast << ',';
		}
		// Generate new Poisson time point;
		x = (rand() + 1.0) / (RAND_MAX + 1.0);
		tLast -= log(x) / rate_;
	}
	last_poisson_time_ = tLast;
	sort(synaptic_driven.begin(), synaptic_driven.end(), compSpike);
}
