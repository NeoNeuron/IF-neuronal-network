//******************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Description: Define class Neuron, structure Spike and NeuronState;
//	Date: 2018-05-30
//******************************
#include<iostream>
#include<cmath>
#include<ctime>
#include<algorithm>
#include "../include/neuron.h"
#include "../include/math_helper.h"
#include "../include/fmath.hpp"
#define exp(x) fmath::expd(x)
using namespace std;

bool compSpike(const Spike &x, const Spike &y) { return x.t < y.t; }

void NeuronSim::GenerateInternalPoisson(bool function, double tmax, bool outSet) {
	double temp, rate, strength;
	if (function) {
		temp = latest_excitatory_poisson_time_;
		rate = excitatory_poisson_rate_;
		strength = feedforward_excitatory_intensity_;
	} else {
		temp = latest_inhibitory_poisson_time_;
		rate = inhibitory_poisson_rate_;
		strength = feedforward_inhibitory_intensity_;
	}
	Spike ADD;
	ADD.t = temp;
	ADD.function = function;
	ADD.s = strength;
	double x, tLast;
	if (rate > 1e-18) {
		if (temp < 1e-18) {
			synaptic_driven_.push_back(ADD);
			if (outSet) cout << ADD.t << '\t';
		}
		tLast = temp;
		while (tLast < tmax) {
			x = (rand() + 1.0) / (RAND_MAX + 1.0);
			tLast -= log(x) / rate;
			// if (tLast > 1e9) cout << "WARNNING: " << tLast << endl;
			ADD.t = tLast;
			synaptic_driven_.push_back(ADD);
			if (outSet) cout << ADD.t << '\t';
		}
		if (function) {
			latest_excitatory_poisson_time_ = tLast;
		} else {
			latest_inhibitory_poisson_time_ = tLast;
		}
	}
	//sort(synaptic_driven_.begin(), synaptic_driven_.end(), compSpike);
}

void NeuronSim::InputExternalPoisson(double tmax, vector<Spike>& x) {
	if (!x.empty()) {
		vector<Spike>::iterator it = x.begin();
		while (it->t < tmax) {
			synaptic_driven_.push_back(*it);
			it = x.erase(it);
			if (x.empty()) break;
		}
	}
	//sort(synaptic_driven_.begin(), synaptic_driven_.end(), compSpike);
}

//void Neuron::SetSynapticStrength(bool function, double S) {
//	if (function)	pyramidal_synaptic_intensity_ = S;
//	else interneuronal_synaptic_intensity_ = S;
//}

void NeuronSim::SetFeedforwardStrength(bool function, double F) {
	if (function)	feedforward_excitatory_intensity_ = F;
	else feedforward_inhibitory_intensity_ = F;
}

void NeuronSim::SetPoissonRate(bool function, double rate) {
	if (function) {
		excitatory_poisson_rate_ = rate;
	} else {
		inhibitory_poisson_rate_ = rate;
	}
}

void NeuronSim::LoadNeuronalState(NeuronalState & data) {
	type_ = data.type;
	//dym_val[v_idx_] = data.membrane_potential_;
	//dym_val[gE_idx_] = data.ge;
	//dym_val[gI_idx_] = data.gi;
	//dym_val[tr_idx_] = data.remaining_refractory_time;
}

void NeuronSim::Reset(double *dym_val) {
	synaptic_driven_.clear();
	spike_train_.clear();
	driven_type_ = false;
	latest_excitatory_poisson_time_ = 0;
	latest_inhibitory_poisson_time_ = 0;
	excitatory_poisson_rate_ = 1e-20;
	inhibitory_poisson_rate_ = 1e-20;
	cycle_ = 0;
	// reset dynamic variables;
	neuron_.ResetDymVal(dym_val);
}

void NeuronSim::OutSpikeTrain(vector<double> & spikes) {
	spikes.clear();
	spikes = spike_train_;
}

void NeuronSim::GetNewSpikes(double t, vector<Spike>& x) {
	Spike add_spike;
	add_spike.s = 0.0;
	add_spike.function = type_;
	x.clear();
	for (vector<double>::reverse_iterator iter = spike_train_.rbegin(); iter != spike_train_.rend(); iter++) {
		if (*iter >= t) {
			add_spike.t = *iter;
			x.push_back(add_spike);
			// cout << *iter << '\t';
		} else break;
		// cout << endl;
	}
}

double NeuronSim::UpdateNeuronalState(double *dym_val, double t, double dt, vector<Spike>& extPoisson, vector<double>& new_spikes) {
	new_spikes.clear();
	double tmax = t + dt;
	if (driven_type_) {
		InputExternalPoisson(tmax, extPoisson);
	} else {
		GenerateInternalPoisson(true, tmax, false);
		//GenerateInternalPoisson(false, tmax, false);
	}
	double t_spike;
	vector<Spike>::iterator s_begin = synaptic_driven_.begin();
	if (s_begin == synaptic_driven_.end() || tmax <= s_begin->t) {
		t_spike = neuron_.DymCore(dym_val, dt);
		cycle_ ++;
		if (t_spike >= 0) new_spikes.push_back(t_spike);
	} else {
		if (t != s_begin->t) {
			t_spike = neuron_.DymCore(dym_val, s_begin->t - t);
			cycle_ ++;
			if (t_spike >= 0) new_spikes.push_back(t_spike);
		}
		for (vector<Spike>::iterator iter = s_begin; iter != synaptic_driven_.end(); iter++) {
			// Update conductance due to the synaptic inputs;
			if (iter -> s) dym_val[ neuron_.GetGEID() ] += iter -> s;
			else dym_val[ neuron_.GetGIID() ] += iter -> s;
			if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
				t_spike = neuron_.DymCore(dym_val, tmax - iter->t);
				cycle_ ++;
				if (t_spike >= 0) new_spikes.push_back(t_spike);
				break;
			} else {
				t_spike = neuron_.DymCore(dym_val, (iter + 1)->t - iter->t);
				cycle_ ++;
				if (t_spike >= 0) new_spikes.push_back(t_spike);
			}
		}
	}
	return dym_val[neuron_.GetVID()];
}

double NeuronSim::CleanUsedInputs(double *dym_val, double *dym_val_new, double tmax) {
	// Update dym_val with dym_val_new;
	for (int i = 0; i < neuron_.GetDymN(); i ++) dym_val[i] = dym_val_new[i];
	// clean old synaptic driven;
	int slen = synaptic_driven_.size();
	if (slen != 0) {
		int i = 0;
		for (; i < slen; i ++) {
			if (synaptic_driven_[i].t >= tmax) break;
		}
		synaptic_driven_.erase(synaptic_driven_.begin(), synaptic_driven_.begin() + i);
	}
	return dym_val[neuron_.GetVID()];
}

void NeuronSim::UpdateConductance(double *dym_val, double t, double dt) {
	double tmax = t + dt;
	if (synaptic_driven_.empty() || tmax <= synaptic_driven_.begin()->t) {
		neuron_.UpdateConductanceOfFiredNeuron(dym_val, dt);
	} else {
		if (t != synaptic_driven_.begin()->t) {
			neuron_.UpdateConductanceOfFiredNeuron(dym_val, synaptic_driven_.begin()->t - t);
		}
		for (vector<Spike>::iterator iter = synaptic_driven_.begin(); iter != synaptic_driven_.end(); iter++) {
			if (iter -> s) dym_val[ neuron_.GetGEID() ] += iter -> s;
			else dym_val[ neuron_.GetGIID() ] += iter -> s;
			if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
				neuron_.UpdateConductanceOfFiredNeuron(dym_val, tmax - iter->t);
				break;
			} else {
				neuron_.UpdateConductanceOfFiredNeuron(dym_val, (iter + 1)->t - iter->t);
			}
		}
	}
}

void NeuronSim::Fire(double *dym_val, double t, double dt) {
	double tmax = t + dt;
	if (synaptic_driven_.empty() || tmax <= synaptic_driven_.begin()->t) {
		neuron_.UpdateConductanceOfFiredNeuron(dym_val, dt);
	} else {
		if (t != synaptic_driven_.begin()->t) {
			neuron_.UpdateConductanceOfFiredNeuron(dym_val, synaptic_driven_.begin()->t - t);
		}
		for (vector<Spike>::iterator iter = synaptic_driven_.begin(); iter != synaptic_driven_.end(); iter++) {
			if (iter -> s) dym_val[ neuron_.GetGEID() ] += iter -> s;
			else dym_val[ neuron_.GetGIID() ] += iter -> s;
			if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
				neuron_.UpdateConductanceOfFiredNeuron(dym_val, tmax - iter->t);
				break;
			} else {
				neuron_.UpdateConductanceOfFiredNeuron(dym_val, (iter + 1)->t - iter->t);
			}
		}
	}
	spike_train_.push_back(t + dt);
	//cout << spike_train_.size() << endl;
	neuron_.ManuallyFire(dym_val);
	int slen = synaptic_driven_.size();
	if (slen != 0) {
		int i = 0;
		for (; i < slen; i ++) {
			if (synaptic_driven_[i].t >= tmax) break;
		}
		synaptic_driven_.erase(synaptic_driven_.begin(), synaptic_driven_.begin() + i);
	}
}

void NeuronSim::Fire(double t, vector<double>& spike_times) {
	for (vector<double>::iterator it = spike_times.begin(); it != spike_times.end(); it ++) {
		spike_train_.push_back(t + *it);
	}
}	

void NeuronSim::InSpike(Spike x) {
	// synaptic_driven_.push_back(x);
	if ((synaptic_driven_.back()).t < x.t) {
		synaptic_driven_.push_back(x);
	} else {
		synaptic_driven_.push_back(x);
		sort(synaptic_driven_.begin(),synaptic_driven_.end(),compSpike);
	}
	// cout << x.t << endl;
}

//double NeuronSim::OutTotalCurrent(double *dym_val) {
//	return - g_m_ * (dym_val[v_idx_] - resting_potential_)
//		- dym_val[gE_idx_] * (dym_val[v_idx_] - excitatory_reversal_potential_)
//		- dym_val[gI_idx_] * (dym_val[v_idx_] - inhibitory_reversal_potential_);
//	return GetDv(dym_val);
//}

//double NeuronSim::OutLeakyCurrent(double *dym_val) {
//	return -g_m_ * (dym_val[v_idx_] - resting_potential_);
//}

//double NeuronSim::OutSynapticCurrent(double *dym_val, bool type) {
//	if (type) {
//		return - dym_val[gE_idx_] * (dym_val[v_idx_] - excitatory_reversal_potential_);
//	} else {
//		return - dym_val[gI_idx_] * (dym_val[v_idx_] - inhibitory_reversal_potential_);
//	}
//}

double NeuronSim::GetConductance(double *dym_val, bool x) {
	if (x) return dym_val[neuron_.GetGEID()];
	else return dym_val[neuron_.GetGIID()];
}

void NeuronSim::Save(NeuronalState & vals) {
	vals.type = type_;
	//	vals.membrane_potential_ = dym_val[v_idx_];
	//	vals.ge = dym_val[gE_idx_];
	//	vals.gi = dym_val[gI_idx_];
	//	vals.remaining_refractory_time = dym_val[tr_idx_];
}

void GenerateExternalPoissonSequence(double rate, double tmax, int seed, vector<double> & list) {
	srand(seed);
	list.push_back(0);
	double x, tLast = 0;
	while (tLast < tmax) {
		x = rand() / (RAND_MAX + 1.0);
		while (x == 0) x = rand() / (RAND_MAX + 1.0);
		tLast -= log(x) / rate;
		list.push_back(tLast);
		//cout << tLast << ',';
	}
	//cout <<endl;
}
