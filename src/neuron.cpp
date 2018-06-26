//******************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Description: Define class Neuron, structure Spike and NeuronState;
//	Date: 2018-05-30
//******************************
#include<iostream>
#include<cmath>
#include<ctime>
#include<vector>
#include<algorithm>
#include "../include/neuron.h"
#include "../include/math_helper.h"
using namespace std;

bool compPoisson(const Spike & x, const Spike & y) {
	return x.t < y.t;
}

void Neuron::GenerateInternalPoisson(bool function, double tmax, bool outSet) {
	double temp, rate;
	if (function == true) {
		temp = latest_excitatory_poisson_time_;
		rate = excitatory_poisson_rate_;
	} else {
		temp = latest_inhibitory_poisson_time_;
		rate = inhibitory_poisson_rate_;
	}
	Spike ADD;
	ADD.t = temp;
	ADD.function = function;
	ADD.mode = true;
	double x, tLast;
	if (rate > 1e-18) {
		if (temp < 1e-18) {
			synaptic_driven_.push_back(ADD);
			if (outSet == true) cout << ADD.t << '\t';
		}
		tLast = temp;
		while (tLast < tmax) {
			x = rand() / (RAND_MAX + 1.0);
			while (x == 0) x = rand() / (RAND_MAX + 1.0);
			tLast -= log(x) / rate;
			// if (tLast > 1e9) cout << "WARNNING: " << tLast << endl;
			ADD.t = tLast;
			synaptic_driven_.push_back(ADD);
			if (outSet == true) cout << ADD.t << '\t';
		}
		if (function == true) {
			latest_excitatory_poisson_time_ = tLast;
		} else {
			latest_inhibitory_poisson_time_ = tLast;
		}
	}
	//sort(synaptic_driven_.begin(), synaptic_driven_.end(), compPoisson);
}

void Neuron::InputExternalPoisson(bool function, double tmax, vector<double>& x) {
	Spike ADD;
	ADD.mode = true;
	ADD.function = function;
	if (x.size() != 0) {
		vector<double>::iterator it = x.begin();
		while (*it < tmax) {
			ADD.t = *it;
			synaptic_driven_.push_back(ADD);
			it = x.erase(it);
			if (x.size() == 0) break;
		}
	}
	//sort(synaptic_driven_.begin(), synaptic_driven_.end(), compPoisson);
}

void Neuron::UpdateG(double *dy_val_, bool mode, bool function) {
	if (mode) {
		if (function) {
			dy_val_[gE_idx_] += feedforward_excitatory_intensity_;
		} else {
			dy_val_[gI_idx_] += feedforward_inhibitory_intensity_;
		}
	} else {
		if (function) {
			dy_val_[gE_idx_] += pyramidal_synaptic_intensity_;
		} else {
			dy_val_[gI_idx_] += interneuronal_synaptic_intensity_;
		}
	}
}

double Neuron::GetDv(double *dy_val_) {
	return - g_m_ * (dy_val_[v_idx_] - resting_potential_)
		- dy_val_[gE_idx_] * (dy_val_[v_idx_] - excitatory_reversal_potential_)
		- dy_val_[gI_idx_] * (dy_val_[v_idx_] - inhibitory_reversal_potential_);
}

double Neuron::UpdatePotential(double *dy_val_, double dt) {
	double exp_E = exp(-0.5 * dt / tau_e_);
	double exp_I = exp(-0.5 * dt / tau_i_);
	// k1 = GetDv(t_n, v_n);
	// k2 = GetDv(t_n+1/2, v_n + k1*dt / 2);
	// k3 = GetDv(t_n+1/2, v_n + k2*dt / 2);
	// k4 = GetDv(t_n+1, v_n + k3*dt);
	// v_n+1 = v_n + dt/6*(k1 + 2*k2 + 2*k3 + k4);
	double v_n = dy_val_[v_idx_];
	double k1, k2, k3, k4;
	k1 = GetDv(dy_val_);
	// Update G:
	dy_val_[gE_idx_] *= exp_E;
	dy_val_[gI_idx_] *= exp_I;
	dy_val_[v_idx_] = v_n + 0.5*k1*dt;
	k2 = GetDv(dy_val_);
	dy_val_[v_idx_] = v_n + 0.5*k2*dt;
	k3 = GetDv(dy_val_);
	// Update G:
	dy_val_[gE_idx_] *= exp_E;
	dy_val_[gI_idx_] *= exp_I;
	dy_val_[v_idx_] = v_n + k3*dt;
	k4 = GetDv(dy_val_);
	// Get v_n+1;
	dy_val_[v_idx_] = v_n + dt / 6 *(k1 + 2 * k2 + 2 * k3 + k4);
	return k1;
}

double Neuron::FindExactSpike(double v1, double dv1, double *dy_val_, double dt) {
	double dv2 = GetDv(dy_val_);
	return cubic_hermite_root(dt, v1, dy_val_[v_idx_], dv1, dv2, threshold_potential_);
}

double Neuron::OffRefractoryPeriod(double *dy_val_, double dt) {
	dy_val_[gE_idx_] *= exp(- remaining_refractory_period_ / tau_e_);
	dy_val_[gI_idx_] *= exp(- remaining_refractory_period_ / tau_i_);
	return UpdatePotential(dy_val_, dt - remaining_refractory_period_);
}

double Neuron::PrimelyUpdateState(double *dy_val_, bool is_fire, bool mode, bool function, double t, double dt, bool temp_switch) {
	//if (temp_switch == false) cout << t << ',';
	double vn = dy_val_[v_idx_];
	// Update conductance;
	if (is_fire) UpdateG(dy_val_, mode, function);
  double dvn, t_spike;
  double return_val = -1; // -1 for no spiking case and non-temp update, otherwise, return_val = t_spike;
	if (remaining_refractory_period_ <= 0) { // neuron is not in the refractory period;
		dvn = UpdatePotential(dy_val_, dt);
		// Check whether fire or not;
		if (dy_val_[v_idx_] >= threshold_potential_) {
			t_spike = FindExactSpike(vn, dvn, dy_val_, dt);
			dy_val_[v_idx_] = resting_potential_;
			//if (temp_switch == false) cout << "mark,";
			//if (temp_switch == false) cout << t_spike << ',';
			//cout << vn << ';';
			if (temp_switch == false) {
				// Add spike time to spike train, update remaining fractory period;
				spike_train_.push_back(t + t_spike);
				//cout << t+t_spike << endl;
				remaining_refractory_period_ = tau_ + t_spike - dt;
			} else return_val = t_spike;
		}	
	} else { // neuron is about to exit the refractory period;
		if (remaining_refractory_period_ < dt) {
			dvn = OffRefractoryPeriod(dy_val_, dt);
			//if (temp_switch == false) cout << remaining_refractory_period_ << ',';
			//if (temp_switch == false) cout << dy_val_[v_idx_] << ',';
			// Check whether fire or not;
			if (dy_val_[v_idx_] >= threshold_potential_) {
				t_spike = FindExactSpike(vn, dvn, dy_val_, dt - remaining_refractory_period_);
				dy_val_[v_idx_] = resting_potential_;
				//cout << t_spike + remaining_refractory_period_ << ',';
				//cout << vn << ';';
				if (temp_switch == false) {
					// Add spike time to spike train, update remaining fractory period;
					spike_train_.push_back(t + t_spike + remaining_refractory_period_);
					//cout << t+t_spike + remaining_refractory_period_ << endl;
					remaining_refractory_period_ = tau_ + t_spike + remaining_refractory_period_;
				} else return_val = t_spike;
			}
		} else { // neuron is in the refractory period;
			dy_val_[gE_idx_] *= exp(- dt / tau_e_);
			dy_val_[gI_idx_] *= exp(- dt / tau_i_);
		}
		remaining_refractory_period_ -= dt;
	}
	//if (temp_switch == false) cout << remaining_refractory_period_ << ',';
	return return_val;
}

void Neuron::UpdateConductanceOfFiredNeuron(double *dy_val_, bool is_fire, bool mode, bool function, double dt) {
	// Update Conductance;
	if (is_fire) UpdateG(dy_val_, mode, function);
	dy_val_[gE_idx_] *= exp(- dt / tau_e_);
	dy_val_[gI_idx_] *= exp(- dt / tau_i_);
}

void Neuron::SetSynapticStrength(bool function, double S) {
	if (function)	pyramidal_synaptic_intensity_ = S;
	else interneuronal_synaptic_intensity_ = S;
}

void Neuron::SetFeedforwardStrength(bool function, double F) {
	if (function)	feedforward_excitatory_intensity_ = F;
	else feedforward_inhibitory_intensity_ = F;
}

void Neuron::SetPoissonRate(bool function, double rate) {
	if (function == true) {
		excitatory_poisson_rate_ = rate;
	} else {
		inhibitory_poisson_rate_ = rate;
	}
}

void Neuron::LoadNeuronalState(NeuronalState & data) {
	type_ = data.type;
	index_ = data.index;
	dy_val_[v_idx_] = data.membrane_potential_;
	dy_val_[gE_idx_] = data.ge;
	dy_val_[gI_idx_] = data.gi;
	remaining_refractory_period_ = data.remaining_refractory_time;
}

void Neuron::Reset() {
	synaptic_driven_.clear();
	spike_train_.clear();
	driven_type_ = false;
	latest_excitatory_poisson_time_ = 0;
	latest_inhibitory_poisson_time_ = 0;
	excitatory_poisson_rate_ = 1e-20;
	inhibitory_poisson_rate_ = 1e-20;
	for (int i = 0; i < 3; i ++) dy_val_[i] = 0.0;
	remaining_refractory_period_ = -1;
}

void Neuron::OutSpikeTrain(vector<double> & spikes) {
	spikes.clear();
	spikes = spike_train_;
}

void Neuron::GetNewSpikes(double t, vector<Spike>& x) {
	Spike add_spike;
	add_spike.mode = false;
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

double Neuron::UpdateNeuronalState(double t, double dt) {
	double tmax = t + dt;
	GenerateInternalPoisson(true, tmax, false);
	//GenerateInternalPoisson(false, tmax, false);
	vector<Spike>::iterator s_begin = synaptic_driven_.begin();
	if (s_begin == synaptic_driven_.end()) {
		PrimelyUpdateState(dy_val_, false, false, false, t, dt, false);
	} else {
		if (tmax <= s_begin->t) {
			PrimelyUpdateState(dy_val_, false, false, false, t, dt, false);
		} else {
			if (t != s_begin->t) {
				PrimelyUpdateState(dy_val_, false, false, false, t, s_begin->t - t, false);
			}
			for (vector<Spike>::iterator iter = s_begin; iter != synaptic_driven_.end(); iter++) {
				if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
					PrimelyUpdateState(dy_val_, true, iter->mode, iter->function, iter->t, tmax - iter->t, false);
					break;
				} else {
					PrimelyUpdateState(dy_val_, true, iter->mode, iter->function, iter->t, (iter + 1)->t - iter->t, false);
				}
			}
		}
	}
	int slen = synaptic_driven_.size();
	if (slen != 0) {
		int i = 0;
		for (; i < slen; i ++) {
			if (synaptic_driven_[i].t >= tmax) break;
		}
		synaptic_driven_.erase(synaptic_driven_.begin(), synaptic_driven_.begin() + i);
	}
	return dy_val_[v_idx_];
}

double Neuron::UpdateNeuronalState(double t, double dt, vector<double> & inPE, vector<double> & inPI) {
	double tmax = t + dt;
	if (driven_type_ == true) {
		InputExternalPoisson(true, tmax, inPE);
		//InputExternalPoisson(false, tmax, inPI);
	} else {
		GenerateInternalPoisson(true, tmax, false);
		//GenerateInternalPoisson(false, tmax, false);
	}
	vector<Spike>::iterator s_begin = synaptic_driven_.begin();
	if (s_begin == synaptic_driven_.end()) {
		PrimelyUpdateState(dy_val_, false, false, false, t, dt, false);
	} else {
		if (tmax <= s_begin->t) {
			PrimelyUpdateState(dy_val_, false, false, false, t, dt, false);
		} else {
			if (t != s_begin->t) {
				PrimelyUpdateState(dy_val_, false, false, false, t, s_begin->t - t, false);
			}
			for (vector<Spike>::iterator iter = s_begin; iter != synaptic_driven_.end(); iter++) {
				if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
					PrimelyUpdateState(dy_val_, true, iter->mode, iter->function, iter->t, tmax - iter->t, false);
					break;
				} else {
					PrimelyUpdateState(dy_val_, true, iter->mode, iter->function, iter->t, (iter + 1)->t - iter->t, false);
				}
			}
		}
	}
	//if (synaptic_driven_.size() != 0) {
	//	vector<Spike>::reverse_iterator rit = synaptic_driven_.rbegin();
	//	for (; rit != synaptic_driven_.rend(); rit ++) {
	//		if (rit -> t < tmax) break;
	//	}
	//	synaptic_driven_.erase(synaptic_driven_.begin(), (rit--).base());
	//}
	int slen = synaptic_driven_.size();
	if (slen != 0) {
		int i = 0;
		for (; i < slen; i ++) {
			if (synaptic_driven_[i].t >= tmax) break;
		}
		synaptic_driven_.erase(synaptic_driven_.begin(), synaptic_driven_.begin() + i);
	}
	return dy_val_[v_idx_];
}

double Neuron::TemporallyUpdateNeuronalState(double t, double dt, vector<double> & inPE, vector<double> & inPI) {
	double tmax = t + dt;
	if (driven_type_ == true) {
		InputExternalPoisson(true, tmax, inPE);
		//InputExternalPoisson(false, tmax, inPI);
	} else {
		GenerateInternalPoisson(true, tmax, false);
		//GenerateInternalPoisson(false, tmax, false);
	}
	// backup dynamical variables:
	double dy_val_bk[3] = {dy_val_[0], dy_val_[1], dy_val_[2]};
	double remaining_refractory_period_bk = remaining_refractory_period_;
	vector<Spike>::iterator s_begin = synaptic_driven_.begin();
	double return_val = -1;
	if (s_begin == synaptic_driven_.end()) {
		return_val = PrimelyUpdateState(dy_val_, false, false, false, t, dt, true);
	} else {
		if (tmax <= s_begin->t) {
			return_val = PrimelyUpdateState(dy_val_, false, false, false, t, dt, true);
		} else {
			if (t != s_begin->t) {
				return_val = PrimelyUpdateState(dy_val_, false, false, false, t, s_begin->t - t, true);
			}
			if (return_val < 0) {
				for (vector<Spike>::iterator iter = s_begin; iter != synaptic_driven_.end(); iter++) {
					if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
						return_val = PrimelyUpdateState(dy_val_, true, iter->mode, iter->function, iter->t, tmax - iter->t, true);
						if (return_val >= 0) return_val += iter->t - t;
						break;
					} else {
						return_val = PrimelyUpdateState(dy_val_, true, iter->mode, iter->function, iter->t, (iter + 1)->t - iter->t, true);
						if (return_val >= 0) {
							return_val += iter->t - t;
							break;
						}
					}
				}
			}
		}
	}
	// restore dynamical variables;
	for (int i = 0; i < 3; i ++) dy_val_[i] = dy_val_bk[i];
	remaining_refractory_period_ = remaining_refractory_period_bk;
	return return_val;
}

void Neuron::Fire(double t, double dt) {
	double tmax = t + dt;
	if (synaptic_driven_.size() == 0) {
		//cout << t << ',';
		UpdateConductanceOfFiredNeuron(dy_val_, false, false, false, dt);
	} else {
		if (tmax <= synaptic_driven_.begin()->t) {
			//cout << t << ',';
			UpdateConductanceOfFiredNeuron(dy_val_, false, false, false, dt);
		} else {
			if (t != synaptic_driven_.begin()->t) {
				//cout << t << ',';
				UpdateConductanceOfFiredNeuron(dy_val_, false, false, false, synaptic_driven_.begin()->t - t);
			}
			for (vector<Spike>::iterator iter = synaptic_driven_.begin(); iter != synaptic_driven_.end(); iter++) {
				if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
					//cout << iter->t << ',';
					UpdateConductanceOfFiredNeuron(dy_val_, true, iter->mode, iter->function, tmax - iter->t);
					break;
				} else {
					//cout << iter->t << ',';
					UpdateConductanceOfFiredNeuron(dy_val_, true, iter->mode, iter->function, (iter + 1)->t - iter->t);
				}
			}
		}
	}
	spike_train_.push_back(t + dt);
	//cout << spike_train_.size() << endl;
	dy_val_[v_idx_] = resting_potential_;
	remaining_refractory_period_ = tau_;
	int slen = synaptic_driven_.size();
	if (slen != 0) {
		int i = 0;
		for (; i < slen; i ++) {
			if (synaptic_driven_[i].t >= tmax) break;
		}
		synaptic_driven_.erase(synaptic_driven_.begin(), synaptic_driven_.begin() + i);
	}
}

void Neuron::InSpike(Spike x) {
	// synaptic_driven_.push_back(x);
	if ((synaptic_driven_.back()).t < x.t) {
		synaptic_driven_.push_back(x);
	} else {
		synaptic_driven_.push_back(x);
		sort(synaptic_driven_.begin(),synaptic_driven_.end(),compPoisson);
	}
	// cout << x.t << endl;
}

double Neuron::OutTotalCurrent() {
	return GetDv(dy_val_);
}

double Neuron::OutLeakyCurrent() {
	return -g_m_ * (dy_val_[v_idx_] - resting_potential_);
}

double Neuron::OutSynapticCurrent(bool type) {
	if (type == true) {
		return - dy_val_[gE_idx_] * (dy_val_[v_idx_] - excitatory_reversal_potential_);
	} else {
		return - dy_val_[gI_idx_] * (dy_val_[v_idx_] - inhibitory_reversal_potential_);
	}
}

double Neuron::GetConductance(bool x) {
	if (x == true) return dy_val_[gE_idx_];
	else return dy_val_[gI_idx_];
}

void Neuron::Save(NeuronalState & vals) {
	vals.type = type_;
	vals.index = index_;
	vals.membrane_potential_ = dy_val_[v_idx_];
	vals.ge = dy_val_[gE_idx_];
	vals.gi = dy_val_[gI_idx_];
	vals.remaining_refractory_time = remaining_refractory_period_;
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
