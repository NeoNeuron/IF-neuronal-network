//******************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Description: Define class Neuron, structure Spike and NeuronState;
//	Date: 2018-04-07
//******************************
#include "../include/neuron.h"
#include<iostream>
#include<cmath>
#include<ctime>
#include<vector>
#include<algorithm>

using namespace std;

bool compPoisson(const Spike & x, const Spike & y) {
	return x.t < y.t;
}

void Neuron::GenerateInternalPoisson(double tmax, bool outSet) {
	double temp, rate;
	temp = latest_excitatory_poisson_time_;
	rate = excitatory_poisson_rate_;
	Spike ADD;
	ADD.t = temp;
	ADD.function = true;
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
		latest_excitatory_poisson_time_ = tLast;
	}
	sort(synaptic_driven_.begin(), synaptic_driven_.end(), compPoisson);
}

void Neuron::InputExternalPoisson(double tmax, vector<double>& x) {
	Spike ADD;
	ADD.mode = true;
	ADD.function = true;
	if (x.size() != 0) {
		vector<double>::iterator it = x.begin();
		while (*it < tmax) {
			ADD.t = *it;
			synaptic_driven_.push_back(ADD);
			it = x.erase(it);
			if (x.size() == 0) break;
		}
	}
	sort(synaptic_driven_.begin(), synaptic_driven_.end(), compPoisson);
}

void Neuron::UpdateConductance(bool mode, bool function, double t, double dt) {
	if (mode == true) {
		if (function == true) {
			excitatory_conductance_tmp_ = excitatory_conductance_tmp_ + feedforward_excitatory_intensity_;
		} else {
			cout << "WARNNING!\n";
		}
	} else {
		if (function == true) {
			excitatory_conductance_tmp_ = excitatory_conductance_tmp_ + pyramidal_synaptic_intensity_;
		} else {
			inhibitory_conductance_tmp_ = inhibitory_conductance_tmp_ + interneuronal_synaptic_intensity_;
		}
	}
}

double Neuron::Alpha(double dt) {
  return g_m_ + excitatory_conductance_tmp_ * exp(-dt / tau_e_) + inhibitory_conductance_tmp_ * exp(-dt / tau_i_);
}

double Neuron::Beta(double dt) {
  return g_m_*resting_potential_ + excitatory_reversal_potential_ * excitatory_conductance_tmp_ * exp(-dt / tau_e_) + inhibitory_reversal_potential_* inhibitory_conductance_tmp_ * exp(-dt / tau_i_);
}

double Neuron::UpdatePotential(double dt) {
	double k1, k2, k3, k4;
	k1 = -Alpha(0.0)* membrane_potential_temp_ + Beta(0.0);
	k2 = -Alpha(dt/2)* (membrane_potential_temp_ + k1*dt / 2) + Beta(dt/2);
	k3 = -Alpha(dt/2)* (membrane_potential_temp_ + k2*dt / 2) + Beta(dt/2);
	k4 = -Alpha(dt)* (membrane_potential_temp_ + k3*dt) + Beta(dt);
	return membrane_potential_temp_ + dt*(k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

double Neuron::CheckFire(double voltage, double dt) {
	if (voltage < threshold_potential_) {
		membrane_potential_temp_ = voltage;
		return -1;
	} else {
		double voltage_1p, voltage_2p;
		voltage_1p = -Alpha(0.0)*membrane_potential_temp_ + Beta(0.0);
		voltage_2p = -Alpha(dt)*voltage + Beta(dt);
		double error = 1, x;
		x = dt*(threshold_potential_ - membrane_potential_temp_) / (voltage - membrane_potential_temp_);
		while (abs(error) > 1e-15) {
			error = ((2 * (membrane_potential_temp_ - voltage) + dt*(voltage_1p + voltage_2p))*pow(x, 3) / pow(dt, 3)
				+ (3 * (voltage - membrane_potential_temp_) - dt*(2 * voltage_1p + voltage_2p))*pow(x, 2) / pow(dt, 2)
				+ voltage_1p*x + membrane_potential_temp_ - threshold_potential_)
				/ ((2 * (membrane_potential_temp_ - voltage) + dt*(voltage_1p + voltage_2p)) * 3 * pow(x, 2) / pow(dt, 3)
					+ (3 * (voltage - membrane_potential_temp_) - dt*(2 * voltage_1p + voltage_2p)) * 2 * x / pow(dt, 2)
					+ voltage_1p);
			x -= error;
		}
		return x;
	}
}

double Neuron::OffRefractoryPeriod(double dt) {
	double voltage, numerator, denominator, alpha1, alpha2, alpha3, beta1, beta2, beta3, t_spike, p1, p2, p3, p4;
	alpha1 = Alpha(0.0);
	alpha2 = Alpha(dt/2);
	alpha3 = Alpha(dt);
	beta1 = Beta(0.0);
	beta2 = Beta(dt/2);
	beta3 = Beta(dt);
	t_spike = remaining_refractory_period_temp_;
	p1 = 1 - 3 * pow(t_spike / dt, 2) + 2 * pow(t_spike / dt, 3);
	p2 = 3 * pow(t_spike / dt, 2) - 2 * pow(t_spike / dt, 3);
	p3 = t_spike - 2 * pow(t_spike, 2) / dt + pow(t_spike, 3) / pow(dt, 2);
	p4 = -pow(t_spike, 2) / dt + pow(t_spike, 3) / pow(dt, 2);
	numerator = resting_potential_ - p3*beta1 - p4*beta3 - (p2 - p4*alpha3)*dt / 6 * (beta1 + 4 * beta2 + beta3
		- (alpha2*beta1 + alpha2*beta2 + alpha3*beta2)*dt + (pow(alpha2, 2)*beta1 / 2 + alpha2*alpha3*beta2 / 2)*pow(dt, 2)
		- pow(alpha2, 2)*alpha3*beta1*pow(dt, 3) / 4);
	denominator = (p1 + p2 - p3*alpha1 - p4*alpha3) + (p2 - p4*alpha3)*dt / 6 * (-(alpha1 + 4 * alpha2 + alpha3)
		+ (alpha1*alpha2 + pow(alpha2, 2) + alpha2*alpha3)*dt- (alpha1*pow(alpha2, 2) / 2 + pow(alpha2, 2)*alpha3 / 2)*pow(dt, 2)
		+ alpha1*pow(alpha2, 2)*alpha3 / 4 * pow(dt, 3));
	membrane_potential_temp_ = numerator / denominator;
	voltage = UpdatePotential(dt);
	return voltage;
}

double Neuron::PrimelyUpdateState(bool is_fire, bool mode, bool function, double t, double dt, bool temp_switch) {
	if (is_fire) {
	  UpdateConductance(mode, function, t, dt);
	}
  double voltage;
  double spike_value = -1; // if negative, no spiking event or non-temp situation; else, spike_value = spiking time node;
	if (temp_switch == false) {
		double spike_time;
		if (remaining_refractory_period_temp_ < 0) {
			voltage = UpdatePotential(dt);
			spike_time = CheckFire(voltage, dt);
			if (spike_time > 0) {
				spike_train_.push_back(spike_time + t);
				remaining_refractory_period_temp_ = tau_ + spike_time - dt;
				membrane_potential_temp_ = resting_potential_;
			}
		} else {
			if (remaining_refractory_period_temp_ < dt) {
				voltage = OffRefractoryPeriod(dt);
				spike_time = CheckFire(voltage, dt - remaining_refractory_period_temp_);
				if (spike_time > 0) {
					spike_train_.push_back(spike_time + t + remaining_refractory_period_temp_);
					remaining_refractory_period_temp_ = tau_ + spike_time - dt + remaining_refractory_period_temp_;
					membrane_potential_temp_ = resting_potential_;
				}
			}
			remaining_refractory_period_temp_ -= dt;
		}
	} else {
		if (remaining_refractory_period_temp_ < 0) {
			voltage = UpdatePotential(dt);
			spike_value = CheckFire(voltage, dt);
		} else {
			if (remaining_refractory_period_temp_ < dt) {
				voltage = OffRefractoryPeriod(dt);
				spike_value = CheckFire(voltage, dt - remaining_refractory_period_temp_);
			}
			remaining_refractory_period_temp_ -= dt;
		}
	}
	excitatory_conductance_tmp_ *= exp(-dt / tau_e_);
	inhibitory_conductance_tmp_ *= exp(-dt / tau_i_);
	return spike_value;
}

void Neuron::UpdateConductanceOfFiredNeuron(bool is_fire, bool mode, bool function, double t, double dt) {
	if (is_fire) {
		UpdateConductance(mode, function, t, dt);
	}
	excitatory_conductance_tmp_ *= exp(-dt / tau_e_);
	inhibitory_conductance_tmp_ *= exp(-dt / tau_i_);
}

void Neuron::SetDrivingType(bool x) {
	driven_type_ = x;
}

void Neuron::SetSynapticStrength(bool function, double S) {
	if (function)	pyramidal_synaptic_intensity_ = S;
	else interneuronal_synaptic_intensity_ = S;
}

void Neuron::SetFeedforwardStrength(double F) {
	feedforward_excitatory_intensity_ = F;
}

void Neuron::SetPoissonRate(double rate) {
	excitatory_poisson_rate_ = rate;
}

void Neuron::LoadNeuronalState(NeuronalState & data) {
	type_ = data.type;
	index_ = data.index;
	membrane_potential_ = data.membrane_potential_;
	excitatory_conductance_ = data.ge;
	inhibitory_conductance_ = data.gi;
	remaining_refractory_period_ = data.remaining_refractory_time;
}

void Neuron::Reset() {
	synaptic_driven_.clear();
	spike_train_.clear();
	driven_type_ = false;
	latest_excitatory_poisson_time_ = 0.0;
	excitatory_poisson_rate_ = 1e-20;
	excitatory_conductance_ = 0.0;
	excitatory_conductance_tmp_ = 0.0;
	inhibitory_conductance_ = 0.0;
	inhibitory_conductance_tmp_ = 0.0;
	membrane_potential_temp_ = 0.0;
	membrane_potential_ = 0.0;
	remaining_refractory_period_temp_ = -1.0;
	remaining_refractory_period_ = -1.0;
}

void Neuron::SetNeuronType(bool x) {
	type_ = x;
}

void Neuron::SetNeuronIndex(int x) {
	index_ = x;
}

double Neuron::GetLastSpike() {
	return spike_train_.back();
}

double Neuron::GetPotential() {
	return membrane_potential_;
}

bool Neuron::GetNeuronalType() {
	return type_;
}

int Neuron::GetNeuronIndex() {
  return index_;
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
	GenerateInternalPoisson(tmax, false);
	excitatory_conductance_tmp_ = excitatory_conductance_;
	inhibitory_conductance_tmp_ = inhibitory_conductance_;
	membrane_potential_temp_ = membrane_potential_;
	remaining_refractory_period_temp_ = remaining_refractory_period_;
	if (synaptic_driven_.size() == 0) {
		PrimelyUpdateState(false, false, false, t, dt, false);
	} else {
		if (tmax < synaptic_driven_.begin()->t) {
			PrimelyUpdateState(false, false, false, t, dt, false);
		} else if (t == synaptic_driven_.begin()->t) {
			for (vector<Spike>::iterator iter = synaptic_driven_.begin(); iter != synaptic_driven_.end(); iter++) {
				if (iter->t >= tmax) break;
				if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
					PrimelyUpdateState(true, iter->mode, iter->function, iter->t, tmax - iter->t, false);
				} else {
					PrimelyUpdateState(true, iter->mode, iter->function, iter->t, (iter + 1)->t - iter->t, false);
				}
			}
		} else if (t< synaptic_driven_.begin()->t && tmax > synaptic_driven_.begin()->t) {
			PrimelyUpdateState(false, false, false, t, synaptic_driven_.begin()->t - t, false);
			for (vector<Spike>::iterator iter = synaptic_driven_.begin(); iter != synaptic_driven_.end(); iter++) {
				if (iter->t >= tmax) break;
				if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
					PrimelyUpdateState(true, iter->mode, iter->function, iter->t, tmax - iter->t, false);
				} else {
					PrimelyUpdateState(true, iter->mode, iter->function, iter->t, (iter + 1)->t - iter->t, false);
				}
			}
		}
	}
	excitatory_conductance_ = excitatory_conductance_tmp_ * exp(-dt / tau_e_);
	inhibitory_conductance_ = inhibitory_conductance_tmp_ * exp(-dt / tau_i_);
	membrane_potential_ = membrane_potential_temp_;
	remaining_refractory_period_ = remaining_refractory_period_temp_;
	if (synaptic_driven_.size() != 0) {
		vector<Spike>::iterator it = synaptic_driven_.begin();
		while (it->t < tmax) {
			it = synaptic_driven_.erase(it);
			if (it == synaptic_driven_.end()) break;
		}
	}
	return membrane_potential_;
}

double Neuron::UpdateNeuronalState(double t, double dt, vector<double> & inPE) {
	double tmax = t + dt;
	if (driven_type_ == true) {
		InputExternalPoisson(tmax, inPE);
	} else {
		GenerateInternalPoisson(tmax, false);
	}
	excitatory_conductance_tmp_ = excitatory_conductance_;
	inhibitory_conductance_tmp_ = inhibitory_conductance_;
	membrane_potential_temp_ = membrane_potential_;
	remaining_refractory_period_temp_ = remaining_refractory_period_;
	if (synaptic_driven_.size() == 0) {
		PrimelyUpdateState(false, false, false, t, dt, false);
	} else {
		if (tmax < synaptic_driven_.begin()->t) {
			PrimelyUpdateState(false, false, false, t, dt, false);
		} else if (t == synaptic_driven_.begin()->t) {
			for (vector<Spike>::iterator iter = synaptic_driven_.begin(); iter != synaptic_driven_.end(); iter++) {
				if (iter->t >= tmax) break;
				if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
					PrimelyUpdateState(true, iter->mode, iter->function, iter->t, tmax - iter->t, false);
				} else {
					PrimelyUpdateState(true, iter->mode, iter->function, iter->t, (iter + 1)->t - iter->t, false);
				}
			}
		} else if (t< synaptic_driven_.begin()->t && tmax > synaptic_driven_.begin()->t) {
			PrimelyUpdateState(false, false, false, t, synaptic_driven_.begin()->t - t, false);
			for (vector<Spike>::iterator iter = synaptic_driven_.begin(); iter != synaptic_driven_.end(); iter++) {
				if (iter->t >= tmax) break;
				if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
					PrimelyUpdateState(true, iter->mode, iter->function, iter->t, tmax - iter->t, false);
				} else {
					PrimelyUpdateState(true, iter->mode, iter->function, iter->t, (iter + 1)->t - iter->t, false);
				}
			}
		}
	}
	excitatory_conductance_ = excitatory_conductance_tmp_ * exp(-dt / tau_e_);
	inhibitory_conductance_ = inhibitory_conductance_tmp_ * exp(-dt / tau_i_);
	membrane_potential_ = membrane_potential_temp_;
	remaining_refractory_period_ = remaining_refractory_period_temp_;
	if (synaptic_driven_.size() != 0) {
		vector<Spike>::iterator it = synaptic_driven_.begin();
		while (it->t < tmax) {
			it = synaptic_driven_.erase(it);
			if (it == synaptic_driven_.end()) break;
		}
	}
	return membrane_potential_;
}

double Neuron::TemporallyUpdateNeuronalState(double t, double dt, vector<double> & inPE) {
	double tmax = t + dt;
	if (driven_type_ == true) {
		InputExternalPoisson(tmax, inPE);
	} else {
		GenerateInternalPoisson(tmax, false);
	}
	excitatory_conductance_tmp_ = excitatory_conductance_;
	inhibitory_conductance_tmp_ = inhibitory_conductance_;
	membrane_potential_temp_ = membrane_potential_;
	remaining_refractory_period_temp_ = remaining_refractory_period_;
	double SET = -1;
	if (synaptic_driven_.size() == 0) {
		SET = PrimelyUpdateState(false, false, false, t, dt, true);
	} else {
		if (tmax < synaptic_driven_.begin()->t) {
			SET = PrimelyUpdateState(false, false, false, t, dt, true);
		} else if (t == synaptic_driven_.begin()->t) {
			for (vector<Spike>::iterator iter = synaptic_driven_.begin(); iter != synaptic_driven_.end(); iter++) {
				if (iter->t >= tmax) break;
				if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
					SET = PrimelyUpdateState(true, iter->mode, iter->function, iter->t, tmax - iter->t, true);
					if (SET > 0) break;
				} else {
					SET = PrimelyUpdateState(true, iter->mode, iter->function, iter->t, (iter + 1)->t - iter->t, true);
					if (SET > 0) break;
				}
			}
		} else if (t< synaptic_driven_.begin()->t && tmax > synaptic_driven_.begin()->t) {
			SET = PrimelyUpdateState(false, false, false, t, synaptic_driven_.begin()->t - t, true);
			if (SET > 0) return SET;
			for (vector<Spike>::iterator iter = synaptic_driven_.begin(); iter != synaptic_driven_.end(); iter++) {
				if (iter->t >= tmax) break;
				if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
					SET = PrimelyUpdateState(true, iter->mode, iter->function, iter->t, tmax - iter->t, true);
					if (SET > 0) break;
				} else {
					SET = PrimelyUpdateState(true, iter->mode, iter->function, iter->t, (iter + 1)->t - iter->t, true);
					if (SET > 0) break;
				}
			}
		}
	}
	return SET;
}

void Neuron::Fire(double t, double dt) {
	double tmax = t + dt;
	excitatory_conductance_tmp_ = excitatory_conductance_;
	inhibitory_conductance_tmp_ = inhibitory_conductance_;
	if (synaptic_driven_.size() == 0) {
		UpdateConductanceOfFiredNeuron(false, false, false, t, dt);
	} else {
		if (tmax < synaptic_driven_.begin()->t || synaptic_driven_.size() == 0) {
			UpdateConductanceOfFiredNeuron(false, false, false, t, dt);
		} else if (t == synaptic_driven_.begin()->t) {
			for (vector<Spike>::iterator iter = synaptic_driven_.begin(); iter != synaptic_driven_.end(); iter++) {
				if (iter->t >= tmax) break;
				if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
					UpdateConductanceOfFiredNeuron(true, iter->mode, iter->function, iter->t, tmax - iter->t);
				} else {
					UpdateConductanceOfFiredNeuron(true, iter->mode, iter->function, iter->t, (iter + 1)->t - iter->t);
				}
			}
		} else if (t< synaptic_driven_.begin()->t && tmax > synaptic_driven_.begin()->t) {
			UpdateConductanceOfFiredNeuron(false, false, false, t, synaptic_driven_.begin()->t - t);
			for (vector<Spike>::iterator iter = synaptic_driven_.begin(); iter != synaptic_driven_.end(); iter++) {
				if (iter->t >= tmax) break;
				if (iter + 1 == synaptic_driven_.end() || (iter + 1)->t >= tmax) {
					UpdateConductanceOfFiredNeuron(true, iter->mode, iter->function, iter->t, tmax - iter->t);
				} else {
					UpdateConductanceOfFiredNeuron(true, iter->mode, iter->function, iter->t, (iter + 1)->t - iter->t);
				}
			}
		}
	}
	spike_train_.push_back(t + dt);
	excitatory_conductance_ = excitatory_conductance_tmp_ * exp(-dt / tau_e_);
	inhibitory_conductance_ = inhibitory_conductance_tmp_ * exp(-dt / tau_i_);
	membrane_potential_ = resting_potential_;
	remaining_refractory_period_ = tau_;
	if (synaptic_driven_.size() != 0) {
		vector<Spike>::iterator it = synaptic_driven_.begin();
		while (it->t < tmax) {
			it = synaptic_driven_.erase(it);
			if (it == synaptic_driven_.end()) break;
		}
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
	return -Alpha(0)*membrane_potential_ + Beta(0);
}

double Neuron::OutLeakyCurrent() {
	return -g_m_ * (membrane_potential_ - resting_potential_);
}

double Neuron::OutSynapticCurrent(bool type) {
	if (type == true) {
		return - excitatory_conductance_ * (membrane_potential_ - excitatory_reversal_potential_);
	} else {
		return - inhibitory_conductance_ * (membrane_potential_ - inhibitory_reversal_potential_);
	}
}

double Neuron::GetConductance(bool x) {
	if (x == true) return excitatory_conductance_;
	else return inhibitory_conductance_;
}

void Neuron::Save(NeuronalState & vals) {
	vals.type = type_;
	vals.index = index_;
	vals.membrane_potential_ = membrane_potential_;
	vals.ge = excitatory_conductance_;
	vals.gi = inhibitory_conductance_;
	vals.remaining_refractory_time = remaining_refractory_period_;
}

void GenerateExternalPoissonSequence(double rate, double tmax, int seed, vector<double> & list) {
	srand(seed);
	list.push_back(0);
	double x, tLast = 0;
	while (tLast < tmax) {
		x = rand() / (RAND_MAX + 1.0);
		tLast -= log(x) / rate;
		list.push_back(tLast);
	}
}
