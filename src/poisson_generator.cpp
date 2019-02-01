#include "../include/poisson_generator.h"
#include <cmath>
#include <cstdlib>

void PoissonGenerator::GenerateNewPoisson( double tmax, queue<Spike>& poisson_driven) {
	Spike new_spike;
	new_spike.type = true;
	new_spike.s = strength_;
	double x, tLast = last_poisson_time_;
	while (tLast < tmax) {
		new_spike.t = tLast;
		poisson_driven.push(new_spike);
		if (output_flag_) {
			outfile_ << setprecision(18) << tLast << ',';
		}
		// Generate new Poisson time point;
		x = (rand() + 1.0) / (RAND_MAX + 1.0);
		tLast -= log(x) / rate_;
	}
	last_poisson_time_ = tLast;
	//sort(poisson_driven.begin(), poisson_driven.end(), compSpike);
}
