//******************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Description: Define class Neuron, structure Spike and NeuronState;
//	Date: 2018-09-27
//******************************
#ifndef _NEURON_H_
#define _NEURON_H_

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "../include/math_helper.h"

using namespace std;

struct Spike {
	bool function; // function of spike: true for excitation(AMPA), false for inhibition(GABA);
	double t; // Exact spiking time;
	double s; // strength of spikes;
};

bool compSpike(const Spike &x, const Spike &y);

struct NeuronalState {
	bool type; // neuronal type: true for excitatory, false for inhibitory;
	int index; // index of neuron in loop lattice;
	double membrane_potential_; // membrane potential;
	double ge; // excitatory conductance;
	double gi; // inhibitory conductance;
	double remaining_refractory_time;
};

// Class Neuron: Based on integrate and fire neuron model;
class LIF_G_Model {
	public:
		// PARAMETERS:
		double tau_e_ = 2.0;	// (ms) time const for excitatory conductance;
		double tau_i_ = 5.0;	// (ms) time const for inhibitory conductance;
		double g_m_ = 5e-2;		// (1/ms) normalized membrane conductance;
		double tau_ = 2.0;		// (ms) refractory Period;
		double resting_potential_ = 0.0;
		double threshold_potential_ = 1.0;
		double excitatory_reversal_potential_ = 14.0 / 3.0;
		double inhibitory_reversal_potential_ = -2.0 / 3.0;

		// excitatory and inhibitory conductance; evolve precisely with the given expression;
		const int dym_n_ = 4;
		const int v_idx_ = 0;
		const int gE_idx_ = 1;
		const int gI_idx_ = 2;
		const int tr_idx_ = 3;
		// index of remaining refractory period time. if negative, remaining refractory period equals to zero;

		// DYNAMICS:

		//	Purely update conductance after single time step dt;
		//	dym_val: dynamical variables;
		//	dt: time step;
		//	return: none;
		void UpdateG(double *dym_val, double dt) {
			dym_val[gE_idx_] *= exp( -dt / tau_e_ );
			dym_val[gI_idx_] *= exp( -dt / tau_i_ );
		}

		// ODE govern the dynamic of IF neuron;
		// dym_val: dynamical variables;
		// return: dV/dt, the derivative of V;
		double GetDv(double *dym_val) {
			return - g_m_ * (dym_val[v_idx_] - resting_potential_)
				- dym_val[gE_idx_] * (dym_val[v_idx_] - excitatory_reversal_potential_)
				- dym_val[gI_idx_] * (dym_val[v_idx_] - inhibitory_reversal_potential_);
		}
		
		//	Update the conductance and membrane potential for t = [t_n, t_n + dt];
		//	Description: 4th-order Runge Kutta integration scheme is applied;
		//	*voltage: pointer of voltage, updated after excecution;
		//	dt: size of time step, unit ms;
		//	return: derivative of membrane potential at t = t(n);
		double DymInplaceRK4(double *dym_val, double dt) {
			double exp_E = exp(-0.5 * dt / tau_e_);
			double exp_I = exp(-0.5 * dt / tau_i_);
			// k1 = GetDv(t_n, v_n);
			// k2 = GetDv(t_n+1/2, v_n + k1*dt / 2);
			// k3 = GetDv(t_n+1/2, v_n + k2*dt / 2);
			// k4 = GetDv(t_n+1, v_n + k3*dt);
			// v_n+1 = v_n + dt/6*(k1 + 2*k2 + 2*k3 + k4);
			double v_n = dym_val[v_idx_];
			double k1, k2, k3, k4;
			k1 = GetDv(dym_val);
			// Update G:
			dym_val[gE_idx_] *= exp_E;
			dym_val[gI_idx_] *= exp_I;
			dym_val[v_idx_] = v_n + 0.5*k1*dt;
			k2 = GetDv(dym_val);
			dym_val[v_idx_] = v_n + 0.5*k2*dt;
			k3 = GetDv(dym_val);
			// Update G:
			dym_val[gE_idx_] *= exp_E;
			dym_val[gI_idx_] *= exp_I;
			dym_val[v_idx_] = v_n + k3*dt;
			k4 = GetDv(dym_val);
			// Get v_n+1;
			dym_val[v_idx_] = v_n + dt / 6 *(k1 + 2 * k2 + 2 * k3 + k4);
			return k1;
		}
};

// class Neuron:
// Implement basic operations for sub-timestep dynamics;
template <class NeuronModel>
class Neuron: public NeuronModel {
	using NeuronModel::tau_e_;	// (ms) time const for excitatory conductance;
	using NeuronModel::tau_i_;	// (ms) time const for inhibitory conductance;
	using NeuronModel::g_m_;		// (1/ms) normalized membrane conductance;
	using NeuronModel::tau_;		// (ms) refractory Period;
	using NeuronModel::resting_potential_;
	using NeuronModel::threshold_potential_;
	using NeuronModel::excitatory_reversal_potential_;
	using NeuronModel::inhibitory_reversal_potential_;
	using NeuronModel::dym_n_;
	using NeuronModel::v_idx_;
	using NeuronModel::gE_idx_;
	using NeuronModel::gI_idx_;
	using NeuronModel::tr_idx_;
	using NeuronModel::UpdateG;
	using NeuronModel::GetDv;
	using NeuronModel::DymInplaceRK4;
	public:
		int GetDymN() { return dym_n_; }
		int GetVID() { return v_idx_; }
		int GetGEID() { return gE_idx_; }
		int GetGIID() { return gI_idx_; }
		int GetTRID() { return tr_idx_; }
		//int Get() {}
		double GetRestingPotential() { return resting_potential_; }
		double GetRefTime() { return tau_; }
		//double Get() {}
		void SetRefTime(double t_ref) { tau_ = t_ref; }
		void ResetDymVal(double* dym_val) {
			dym_val[v_idx_] = 0.0;
			dym_val[gE_idx_] = 0.0;
			dym_val[gI_idx_] = 0.0;
			dym_val[tr_idx_] = -1;
		}
		void ManuallyFire(double* dym_val) {
			dym_val[v_idx_] = 0.0;
			dym_val[tr_idx_] = tau_;
		}
		//	Core operation for updating neuronal state within single timing step dt;
		//	Description: operation to update neuronal state in primary level, including updating conductances, membrane potential and checking spiking events; 
		//	dym_val: array of dynamic variables;
		//  dt: size of time step, unit ms;
		//	return: -1 for no spiking events; otherwise, return relative spiking time respect to the begining of the time step;
		double DymCore(double *dym_val, double dt) {
			double vn = dym_val[v_idx_];
			// Update conductance;
			double dvn, dv_new;
			double t_spike = -1; // spike time within dt;
			if (dym_val[tr_idx_] <= 0) { // neuron is not in the refractory period;
				dvn = DymInplaceRK4(dym_val, dt);
				// Check whether fire or not;
				if (dym_val[v_idx_] >= threshold_potential_) {
					dv_new = GetDv(dym_val);
					t_spike = cubic_hermite_root(dt, vn, dym_val[v_idx_], dvn, dv_new, threshold_potential_);
					dym_val[v_idx_] = resting_potential_;
					// update remaining fractory period
					dym_val[tr_idx_] = tau_ + t_spike - dt;
				}	
			} else { // neuron is about to exit the refractory period;
				if (dym_val[tr_idx_] < dt) {
					UpdateG(dym_val, dym_val[tr_idx_]);
					dvn = DymInplaceRK4(dym_val, dt - dym_val[tr_idx_]);
					// Check whether fire or not;
					if (dym_val[v_idx_] >= threshold_potential_) {
						dv_new = GetDv(dym_val);
						t_spike = cubic_hermite_root(dt - dym_val[tr_idx_], vn, dym_val[v_idx_], dvn, dv_new, threshold_potential_);
						dym_val[v_idx_] = resting_potential_;
						// update remaining fractory period
						t_spike += dym_val[tr_idx_];
						dym_val[tr_idx_] = tau_ + t_spike;
					}
				} else { // neuron is in the refractory period;
					UpdateG(dym_val, dt);
				}
				dym_val[tr_idx_] -= dt;
			}
			return t_spike;
		}
		//	Update conductance of fired neuron within single time step dt; it has the same hierachy level as the PrimelyUpdateState(double*, bool, Spike, double, bool);
		//	Description: operation to update neuronal state in primary level, ONE synaptic input most which arrives at the begining of time step;
		//  is_fire: true for the presence of a synaptic input at the begining of time step; false for not;
		//  mode: mode for the synaptic input;
		//  function: function for synaptic input;
		//  dt: size of time step, unit millisecond;
		//	return: none;
	void UpdateConductanceOfFiredNeuron(double *dym_val, double dt) {
		UpdateG(dym_val, dt);
	}

};

typedef Neuron<LIF_G_Model> LIF_G;

// Class Neuron: Based on integrate and fire neuron model;
class NeuronSim {
	private:
		// Neuron;
		LIF_G neuron_;

		// PARAMETERS:
		bool type_;						// neuronal type: true for excitatory, false for inhibitory;
		double feedforward_excitatory_intensity_; // for internal Poisson generators (default: 5e-3);
		double feedforward_inhibitory_intensity_; // for internal Poisson generators (default: 5e-3);

		// DATA:

		size_t cycle_;	// number of cycle that neuron processed;
		vector<double> spike_train_; // Exact time nodes that neuron fires.
		bool driven_type_; // True for external Poisson driven, false for internal Poisson driven;
		double excitatory_poisson_rate_;
		double inhibitory_poisson_rate_;
		// Synaptic input received by neuron, including feedforward and interneuronal spikes;
		vector<Spike> synaptic_driven_;  
		// record last Poisson spike time generated by Poisson generator;
		double latest_excitatory_poisson_time_; 
		double latest_inhibitory_poisson_time_; 

		// FUNCTIONS:

		// Generate Poisson sequence within each time step; autosort after generatation if synaptic delay is nonzero;
		// function: function of Poisson spike, true for exc, false for inh;
		// tmax: maximum time of Poisson sequence;
		// outSet: whether print spike times of each spiking events, true for print, false for not;
		// return: none;
		void GenerateInternalPoisson(bool function, double tmax, bool outSet);

		// Input external Poisson sequence within each time step, autosort after generatation if synaptic delay is nonzero;
		// function: function of Poisson spike, true for exc, false for inh;
		// tmax: maximum time of Poisson sequence;
		// x: container of external inputing spikes;
		// return: none;
		void InputExternalPoisson(double tmax, vector<Spike> & x);

		public:
		// Initialization of parameters in Neuron;
		NeuronSim(double *&dym_val) {
			type_ = true;
			feedforward_excitatory_intensity_ = 5e-3;
			feedforward_inhibitory_intensity_ = 5e-3;
			driven_type_ = false;
			excitatory_poisson_rate_ = 1e-18;
			inhibitory_poisson_rate_ = 1e-18;
			latest_excitatory_poisson_time_ = 0.0;
			latest_inhibitory_poisson_time_ = 0.0;
			cycle_ = 0;
			dym_val = new double[ neuron_.GetDymN() ];
			neuron_.ResetDymVal(dym_val);
		}

		// INPUTS:
		// Set neuronal type: true for excitatory; false for inhibitory;
		void SetNeuronType(bool x) { type_ = x; }

		// Set refractory period:
		void SetRef(double t_ref) { neuron_.SetRefTime(t_ref); }

		//	Set driving type: true for external, false for internal;
		void SetDrivingType(bool x) { driven_type_ = x; }

		// Set neuronal interacting strength;
		//void SetSynapticStrength(bool function, double S);
		void SetFeedforwardStrength(bool function, double F);

		//	Set Poisson Rate: homogeneous Poisson driving rate of internal driving type;
		//	function: type of Poisson drive, true for excitatory, false for inhibitory;
		void SetPoissonRate(bool function, double rate);

		// Define a 'neuron_file' type to store neuronal condition;
		// A ROW VECTOR:
		//	0: neuronal type;
		//	1: neuronal index;
		//	2: membrane potential;
		//	3: excitatory conductivity;
		//	4: inhibitory conductivity;
		//	5: remaining refractory period;
		void LoadNeuronalState(NeuronalState & data);

		//	Input synaptic inputs, either feedforward or interneuronal ones, autosort after insertion;
		void InSpike(Spike x);

		// Reset neuron into the condition at zero time point;
		void Reset(double *dym_val);

		// DYNAMICS:

		// 	Update neuronal state:
		//	Description: update neuron within single time step, including its membrane potential, conductances and counter of refractory period;
		//	DOUBLE t: time point of the begining of the time step;
		//	DOUBLE dt: size of time step;
		//	VECTOR<DOUBLE> inPE: external excitatory Poisson sequence;
		//	VECTOR<DOUBLE> inPI: external inhibitory Poisson sequence;
		//	Return: membrane potential at t = t + dt;
		//double UpdateNeuronalState(double *dym_val, double t, double dt);
		double UpdateNeuronalState(double *dym_val, double t, double dt, vector<Spike> & extPoisson, vector<double>& new_spikes);
		double CleanUsedInputs(double *dym_val, double * dym_val_new, double tmax);

		void UpdateConductance(double *dym_val, double t, double dt);
		//	Fire: update neuronal state for neurons which fire at t = t + dt;
		void Fire(double *dym_val, double t, double dt);
		void Fire(double dt, vector<double>& spike_times);

		// OUTPUTS:

		// Print cycle_:
		size_t GetCycle() {
			cout << cycle_;
			return cycle_;
		}
		//	Get last spike: return the time point of latest spiking events;
		double GetLastSpike() { return spike_train_.back(); }

		//	Get potential: return the current value of membrane potential;
		double GetPotential(double *dym_val) { return dym_val[neuron_.GetVID()]; }

		//	Get neuronal type: true for excitatory, false for inhibitory;
		bool GetNeuronalType() { return type_; }

		//	Output spike train
		void OutSpikeTrain(vector<double> & spikes);

		//  Output Spikes after t;
		//  the interacting strength of Spikes are set as default(0.0);
		void GetNewSpikes(double t, vector<Spike> &x);

		// Total membrane current;
		double OutTotalCurrent(double *dym_val);

		// Leaky current;
		double OutLeakyCurrent(double *dym_val);

		// Excitatory or inhibitory membrane current;
		double OutSynapticCurrent(double *dym_val, bool type);

		// True return excitatory conductance, false return inhibitory conductance;
		double GetConductance(double *dym_val, bool x);

		//Save corrent neuronal States:
		//Define a 'neuronFile' type to store neuronal condition;
		//A ROW VECTOR:
		//	0: neuronal type;
		//	1: neuronal index;
		//	2: membrane potential;
		//	3: excitatory conductivity;
		//	4: inhibitory conductivity;
		//	5: remaining refractory period;
		void Save(NeuronalState & vals);
};

//	external Poisson generator:
//	rate: mean Poisson firing rate;
//	tmax: maximum timel'
//	seed: random seed;
//	list: memory storage for Poisson squence;
void GenerateExternalPoissonSequence(double rate, double tmax, int seed, vector<double> & list);

#endif 	// _NEURON_H_
