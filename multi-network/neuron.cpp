//******************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Description: Define class Neuron, structure Spike and NeuronState;
//	Date: 2017-03-08 11:09:04
//******************************
#include "neuron.h"
#include<iostream>
#include<cmath>
#include<ctime>
#include<vector>
#include<algorithm>

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
			tLast -= log(x) / rate;
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
	sort(synaptic_driven_.begin(), synaptic_driven_.end(), compPoisson);		
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
	sort(synaptic_driven_.begin(), synaptic_driven_.end(), compPoisson);
}

void Neuron::UpdateExcitatoryConductance(bool mode, bool function, double t, double dt) {
	if (mode == true) {
		if (function == true) {
			excitatory_conductance_1_ = excitatory_conductance_1_ + feedforward_excitatory_intensity_;
			excitatory_conductance_2_ = excitatory_conductance_1_ * exp(-dt / tau_e_ / 2);
			excitatory_conductance_3_ = excitatory_conductance_1_ * exp(-dt / tau_e_);
		} else {
			excitatory_conductance_1_ = excitatory_conductance_1_;
			excitatory_conductance_2_ = excitatory_conductance_1_ * exp(-dt / tau_e_ / 2);
			excitatory_conductance_3_ = excitatory_conductance_1_ * exp(-dt / tau_e_);
		}
	} else {
		if (function == true) {
			excitatory_conductance_1_ = excitatory_conductance_1_ + pyramidal_synaptic_intensity_;
			excitatory_conductance_2_ = excitatory_conductance_1_ * exp(-dt / tau_e_ / 2);
			excitatory_conductance_3_ = excitatory_conductance_1_ * exp(-dt / tau_e_);
		} else {
			excitatory_conductance_1_ = excitatory_conductance_1_;
			excitatory_conductance_2_ = excitatory_conductance_1_ * exp(-dt / tau_e_ / 2);
			excitatory_conductance_3_ = excitatory_conductance_1_ * exp(-dt / tau_e_);
		}
	}
}

void Neuron::UpdateInhibitoryConductance(bool mode, bool function, double t, double dt) {
	if (mode == true) {
		if (function == false) {
			inhibitory_conductance_1_ = inhibitory_conductance_1_ + feedforward_inhibitory_intensity_;
			inhibitory_conductance_2_ = inhibitory_conductance_1_ * exp(-dt / tau_i_ / 2);
			inhibitory_conductance_3_ = inhibitory_conductance_1_ * exp(-dt / tau_i_);
		} else {
			inhibitory_conductance_1_ = inhibitory_conductance_1_;
			inhibitory_conductance_2_ = inhibitory_conductance_1_ * exp(-dt / tau_i_ / 2);
			inhibitory_conductance_3_ = inhibitory_conductance_1_ * exp(-dt / tau_i_);
		}
	} else {
		if (function == false) {
			inhibitory_conductance_1_ = inhibitory_conductance_1_ + interneuronal_synaptic_intensity_;
			inhibitory_conductance_2_ = inhibitory_conductance_1_ * exp(-dt / tau_i_ / 2);
			inhibitory_conductance_3_ = inhibitory_conductance_1_ * exp(-dt / tau_i_);
		} else {
			inhibitory_conductance_1_ = inhibitory_conductance_1_;
			inhibitory_conductance_2_ = inhibitory_conductance_1_ * exp(-dt / tau_i_ / 2);
			inhibitory_conductance_3_ = inhibitory_conductance_1_ * exp(-dt / tau_i_);
		}		
	}
}

double Neuron::Alpha(int option) {
  double value;
  switch (option) {
  case 1:
    value = g_m_ + excitatory_conductance_1_ + inhibitory_conductance_1_;
    break;
  case 2:
    value = g_m_ + excitatory_conductance_2_ + inhibitory_conductance_2_;
    break;
  case 3:
    value = g_m_ + excitatory_conductance_3_ + inhibitory_conductance_3_;
    break;
  default:
    value = 0;
  }
  return value;
}

double Neuron::Beta(int option) {
  double value;
  switch (option) {
  case 1:
    value = g_m_*resting_potential_ + excitatory_reversal_potential_*excitatory_conductance_1_ + inhibitory_reversal_potential_*inhibitory_conductance_1_;
    break;
  case 2:
    value = g_m_*resting_potential_ + excitatory_reversal_potential_*excitatory_conductance_2_ + inhibitory_reversal_potential_*inhibitory_conductance_2_;
    break;
  case 3:
    value = g_m_*resting_potential_ + excitatory_reversal_potential_*excitatory_conductance_3_ + inhibitory_reversal_potential_*inhibitory_conductance_3_;
    break;
  default:
    value = 0;
  }
  return value;
}

double Neuron::UpdatePotential(double dt) {
	double k1, k2, k3, k4;
	k1 = -Alpha(1)* membrane_potential_temp_ + Beta(1);
	k2 = -Alpha(2)* (membrane_potential_temp_ + k1*dt / 2) + Beta(2);
	k3 = -Alpha(2)* (membrane_potential_temp_ + k2*dt / 2) + Beta(2);
	k4 = -Alpha(3)* (membrane_potential_temp_ + k3*dt) + Beta(3);
	return membrane_potential_temp_ + dt*(k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

void Neuron::CheckFire(double voltage, double t, double dt) {
	if (voltage < threshold_potential_) { membrane_potential_temp_ = voltage; }
	else {
		double voltage_1p, voltage_2p;
		voltage_1p = -Alpha(1)*membrane_potential_temp_ + Beta(1);
		voltage_2p = -Alpha(3)*voltage + Beta(3);
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
		spike_train_.push_back(x + t);
		remaining_refractory_period_temp_ = tau_ + x - dt;
		membrane_potential_temp_ = resting_potential_;
	}
}

double Neuron::TempCheckFire(double voltage, double t, double dt) {
	if (voltage < threshold_potential_) {  
		membrane_potential_temp_ = voltage;
		return -1;
	} else {
		double voltage_1p, voltage_2p;
		voltage_1p = -Alpha(1)*membrane_potential_temp_ + Beta(1);
		voltage_2p = -Alpha(3)*voltage + Beta(3);
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
		return t + x;
	}
}

double Neuron::OffRefractoryPeriod(double dt) {
	double voltage, numerator, denominator, alpha1, alpha2, alpha3, beta1, beta2, beta3, t_spike, p1, p2, p3, p4;
	alpha1 = Alpha(1);
	alpha2 = Alpha(2);
	alpha3 = Alpha(3);
	beta1 = Beta(1);
	beta2 = Beta(2);
	beta3 = Beta(3);
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
	if (is_fire == false) {
		UpdateExcitatoryConductance(mode, false, t, dt);
		UpdateInhibitoryConductance(mode, true, t, dt);
	}	else {
	  UpdateExcitatoryConductance(mode, function, t, dt);
		UpdateInhibitoryConductance(mode, function, t, dt);
	}
  double voltage;
  double spike_value = -1; // if negative, no spiking event or non-temp situation; else, spike_value = spiking time node;
	if (temp_switch == false) {
		if (remaining_refractory_period_temp_ < 0) {
			voltage = UpdatePotential(dt);
			CheckFire(voltage, t, dt);
		} else {
			if (remaining_refractory_period_temp_ < dt) {
				voltage = OffRefractoryPeriod(dt);
				CheckFire(voltage, t + remaining_refractory_period_temp_, dt - remaining_refractory_period_temp_);
			}
			remaining_refractory_period_temp_ -= dt;
		}
	} else {
		if (remaining_refractory_period_temp_ < 0) {
			voltage = UpdatePotential(dt);
			spike_value = TempCheckFire(voltage, t, dt);
		} else {
			if (remaining_refractory_period_temp_ < dt) {
				voltage = OffRefractoryPeriod(dt);
				spike_value = TempCheckFire(voltage, t + remaining_refractory_period_temp_, dt - remaining_refractory_period_temp_); 
			}
			remaining_refractory_period_temp_ -= dt;
		}
	}
	excitatory_conductance_1_ = excitatory_conductance_3_;
	inhibitory_conductance_1_ = inhibitory_conductance_3_;	
	return spike_value;
}

void Neuron::UpdateConductanceOfFiredNeuron(bool is_fire, bool mode, bool function, double t, double dt) {
	if (is_fire == false) {
		UpdateExcitatoryConductance(mode, false, t, dt);
		UpdateInhibitoryConductance(mode, true, t, dt);
	} else {
		UpdateExcitatoryConductance(mode, function, t, dt);
		UpdateInhibitoryConductance(mode, function, t, dt);
	}
	excitatory_conductance_1_ = excitatory_conductance_3_;
	inhibitory_conductance_1_ = inhibitory_conductance_3_;	
}

void Neuron::SetDrivingType(bool x) {
	driven_type_ = x;
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
	membrane_potential_ = data.membrane_potential_;
	excitatory_conductance_ = data.ge;
	inhibitory_conductance_ = data.gi;
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
	excitatory_conductance_ = 0;
	excitatory_conductance_1_ = 0;
	excitatory_conductance_2_ = 0;
	excitatory_conductance_3_ = 0;
	inhibitory_conductance_ = 0;
	inhibitory_conductance_1_ = 0;
	inhibitory_conductance_2_ = 0;
	inhibitory_conductance_3_ = 0;
	membrane_potential_temp_ = 0;
	membrane_potential_ = 0;
	remaining_refractory_period_temp_ = -1;
	remaining_refractory_period_ = -1;
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

void Neuron::OutputSpikeTrain(vector<double>& x) {
	for (vector<double>::iterator iter = spike_train_.begin(); iter != spike_train_.end(); iter++) {
		x.push_back(*iter);
	}
}

void Neuron::OutputNewSpikes(double t, vector<Spike>& x) {
	Spike add_spike;
	add_spike.mode = false;
	add_spike.function = type_;
	for (vector<double>::reverse_iterator iter = spike_train_.rbegin(); iter != spike_train_.rend(); iter++) {
		if (*iter >= t) {
			add_spike.t = *iter;
			x.push_back(add_spike);
		} else break;
	}
}

double Neuron::UpdateNeuronalState(double t, double dt, vector<double> & inPE, vector<double> & inPI) {
	double tmax = t + dt;
	if (driven_type_ == true) {
		InputExternalPoisson(true, tmax, inPE);
		InputExternalPoisson(false, tmax, inPI);
	} else {
		GenerateInternalPoisson(true, tmax, false);
		GenerateInternalPoisson(false, tmax, false);
	}
	excitatory_conductance_1_ = excitatory_conductance_;
	inhibitory_conductance_1_ = inhibitory_conductance_;
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
	excitatory_conductance_ = excitatory_conductance_3_;
	inhibitory_conductance_ = inhibitory_conductance_3_;
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

double Neuron::TemporallyUpdateNeuronalState(double t, double dt, vector<double> & inPE, vector<double> & inPI) {
	double tmax = t + dt;
	if (driven_type_ == true) {
		InputExternalPoisson(true, tmax, inPE);
		InputExternalPoisson(false, tmax, inPI);
	} else {
		GenerateInternalPoisson(true, tmax, false);
		GenerateInternalPoisson(false, tmax, false);
	}
	excitatory_conductance_1_ = excitatory_conductance_;
	inhibitory_conductance_1_ = inhibitory_conductance_;
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
	excitatory_conductance_1_ = excitatory_conductance_;
	inhibitory_conductance_1_ = inhibitory_conductance_;
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
	excitatory_conductance_ = excitatory_conductance_3_;
	inhibitory_conductance_ = inhibitory_conductance_3_;
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

void Neuron::InputSpike(Spike x) {
	synaptic_driven_.push_back(x);	
}

void Neuron::SetFeedforwardConductance(bool function, double F) {
	if (function == true)	{
		feedforward_excitatory_intensity_ = F;
	} else {
		feedforward_inhibitory_intensity_ = F;
	}
}

double Neuron::OutputIonCurrent() {
	return -Alpha(3)*membrane_potential_ + Beta(3);
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
