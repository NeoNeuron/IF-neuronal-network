//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2018-01-30
//	Description: source file of lfp.h
//***************
#include "../include/lfp.h"
#include "../include/io.h"
#include "../include/common_header.h"
#include "../include/neuron.h"
using namespace std;

bool comp(const int x, const int y) {
	return x < y;
}

inline double L2(vector<double>& a, vector<double>& b) {
	return sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1]));
}

void Sample(vector<int> & origin_vector, vector<int> & sample_vector, int num) {
	sample_vector.clear();
	if (num == origin_vector.size() || num == 0) {
		sample_vector = origin_vector;
		return;
	} else {
		srand(time(0));
		int ind;
		sample_vector.resize(num);
		vector<int> vector_copy = origin_vector;
		for (int i = 0; i < num; i++) {
			if (vector_copy.size() == 1) ind = 0;
			else ind = rand() % (vector_copy.size() - 1);
			sample_vector[i] = vector_copy[ind];
			vector_copy.erase(vector_copy.begin() + ind, vector_copy.begin() + ind + 1);
		}
		return;
	}
}

int KeySelect(string & key, vector<neuron_type> & types, vector<int> & indices) {
	// select keys:
	indices.clear();
	if (key == "all") {
		indices.resize(types.size());
		for (int i = 0; i < types.size(); i ++) {
			indices[i] = types[i].index;
		}
	} else if (key == "exc") {
		for (vector<neuron_type>::iterator it = types.begin(); it != types.end(); it ++) {
			if (it->type == true) indices.push_back(it->index);
		}
	} else if (key == "inh") {
		for (vector<neuron_type>::iterator it = types.begin(); it != types.end(); it ++) {
			if (it->type == false) indices.push_back(it->index);
		}
	}
	return indices.size();
}

void CalculateSpatialWeight(map<string, string> & m_config, vector<double> & spatial_weights) {
	bool is_spatial;
	istringstream(m_config["IsSpatial"]) >> boolalpha >> is_spatial;
	spatial_weights.clear();
	int neuron_num = atoi(m_config["NeuronNumber"].c_str());
	spatial_weights.resize(neuron_num, 1.0);
	if (is_spatial) {
		// Calculate the spatial weights of selected neurons:
		vector<double> electrode_pos = {atof(m_config["PosX"].c_str()), atof(m_config["PosY"].c_str())};
		// read the coordinate file;
		vector<vector<double> > coordinates;
		Read2D(m_config["CoorPath"], coordinates);
		double distance;
		int decay_order = atoi(m_config["DecayOrder"].c_str());
		if (decay_order == 1) {
			for (int i = 0; i < neuron_num; i ++) {
				distance = L2(electrode_pos, coordinates[i]);
				spatial_weights[i] = 1.0 / distance;
			}
		} else if (decay_order == 2) {
			for (int i = 0; i < neuron_num; i ++) {
				distance = L2(electrode_pos, coordinates[i]);
				spatial_weights[i] = 1.0 / distance / distance;
			}
		} else {
			cout << "Not proper decay order\n";
		}
	}
	return;
}

void CalculateLFP(string dir, vector<double>& lfp, vector<int>& neuron_list, string LFP_type, vector<double>& spatial_weights, double* t_range, double sampling_dt) {
	// preliminary parameters;
	double sampling_rate = 1 / sampling_dt; // Unit ms;

	// Preparing time series;
	size_t t_begin = t_range[0] * sampling_rate; // not included
	size_t t_end = t_range[1] * sampling_rate; // included
	size_t size_of_lfp = t_end - t_begin;
	lfp.clear();
	lfp.resize(size_of_lfp);

	// Load neuron model;
	LIF_G_Model nrn_model;
	double g_m = nrn_model.g_m_;
	double V_rest = nrn_model.resting_potential_;
	double V_e = nrn_model.excitatory_reversal_potential_;
	double V_i = nrn_model.inhibitory_reversal_potential_;
	int v_id = nrn_model.v_idx_;
	int ge_id = nrn_model.gE_idx_;
	int gi_id = nrn_model.gI_idx_;

	ifstream V_file, GE_file, GI_file;
	size_t shape[2];
	V_file.open(dir + "/V.bin", ios::binary);
	V_file.read((char*)&shape, 2*sizeof(size_t));
	vector<double> buffer_vec(shape[0]*shape[1]);
	// Preparing diff_list:
	vector<int> diff_list(neuron_list.size() + 1);
	// Sort neuron list:
	if (neuron_list.size() == 1) {
		diff_list[0] = neuron_list[0];
		diff_list[1] = shape[1] - neuron_list[0] - 1;
	} else {
		sort(neuron_list.begin(), neuron_list.end(), comp);
		for (int i = 0; i < diff_list.size(); i ++) {
			if (i == 0) diff_list[i] = neuron_list[0];
			else if (i == diff_list.size() - 1) diff_list[i] = shape[1] - neuron_list[i - 1] - 1;
			else diff_list[i] = neuron_list[i] - neuron_list[i - 1] - 1;
		}
	}

	// Classify the LFP type and load other dym data file;
	if (LFP_type == "tot") {
		GE_file.open(dir + "/GE.bin", ios::binary);
		GI_file.open(dir + "/GI.bin", ios::binary);
		GE_file.read((char*)buffer_vec.data(), 2*sizeof(size_t));
		GI_file.read((char*)buffer_vec.data(), 2*sizeof(size_t));

		// For t = [0, t_begin];
		V_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
		GE_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
		GI_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
		// For t = (t_begin, t_end]
		double dym_val[3] = {0.0, 0.0, 0.0};
		double temp_lfp;
		for (size_t i = t_begin; i < t_end; i++) {
			temp_lfp = 0;
			if (neuron_list.size() == shape[1]) {
				for (int j = 0; j < neuron_list.size(); j ++) {
					V_file.read((char*)&dym_val[0], sizeof(double));
					GE_file.read((char*)&dym_val[1], sizeof(double));
					GI_file.read((char*)&dym_val[2], sizeof(double));
					temp_lfp += (- g_m * (dym_val[v_id] - V_rest) - dym_val[ge_id] * (dym_val[v_id] - V_e) - dym_val[gi_id] * (dym_val[v_id] - V_i)) * spatial_weights[j];
				}
				lfp[i - t_begin] = temp_lfp / neuron_list.size();
			} else {
				for (int j = 0; j < diff_list.size() - 1; j ++) {
					V_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
					GE_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
					GI_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
					V_file.read((char*)&dym_val[0], sizeof(double));
					GE_file.read((char*)&dym_val[1], sizeof(double));
					GI_file.read((char*)&dym_val[2], sizeof(double));
					temp_lfp += (- g_m * (dym_val[v_id] - V_rest) - dym_val[ge_id] * (dym_val[v_id] - V_e) - dym_val[gi_id] * (dym_val[v_id] - V_i)) * spatial_weights[j];
				}
				V_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
				GE_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
				GI_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
				lfp[i - t_begin] = temp_lfp / neuron_list.size();
			}
		}
		GE_file.close();
		GI_file.close();
	}	else if (LFP_type == "lek") {
		// For t = [0, t_begin];
		V_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
		// For t = (t_begin, t_end]
		double dym_val[3] = {0.0, 0.0, 0.0};
		double temp_lfp;
		for (size_t i = t_begin; i < t_end; i++) {
			temp_lfp = 0;
			if (neuron_list.size() == shape[1]) {
				for (int j = 0; j < neuron_list.size(); j ++) {
					V_file.read((char*)&dym_val[0], sizeof(double));
					temp_lfp += (- g_m * (dym_val[v_id] - V_rest)) * spatial_weights[j];
				}
				lfp[i - t_begin] = temp_lfp / neuron_list.size();
			} else {
				for (int j = 0; j < diff_list.size() - 1; j ++) {
					V_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
					V_file.read((char*)&dym_val[0], sizeof(double));
					temp_lfp += (- g_m * (dym_val[v_id] - V_rest)) * spatial_weights[j];
				}
				V_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
				lfp[i - t_begin] = temp_lfp / neuron_list.size();
			}
		}
	} else if (LFP_type == "exi") {
		GE_file.open(dir + "/GE.bin", ios::binary);
		GE_file.read((char*)buffer_vec.data(), 2*sizeof(size_t));

		// For t = [0, t_begin];
		V_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
		GE_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
		// For t = (t_begin, t_end]
		double dym_val[3] = {0.0, 0.0, 0.0};
		double temp_lfp;
		for (size_t i = t_begin; i < t_end; i++) {
			temp_lfp = 0;
			if (neuron_list.size() == shape[1]) {
				for (int j = 0; j < neuron_list.size(); j ++) {
					V_file.read((char*)&dym_val[0], sizeof(double));
					GE_file.read((char*)&dym_val[1], sizeof(double));
					temp_lfp += (- dym_val[ge_id] * (dym_val[v_id] - V_e)) * spatial_weights[j];
				}
				lfp[i - t_begin] = temp_lfp / neuron_list.size();
			} else {
				for (int j = 0; j < diff_list.size() - 1; j ++) {
					V_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
					GE_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
					V_file.read((char*)&dym_val[0], sizeof(double));
					GE_file.read((char*)&dym_val[1], sizeof(double));
					temp_lfp += (- dym_val[ge_id] * (dym_val[v_id] - V_e)) * spatial_weights[j];
				}
				V_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
				GE_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
				lfp[i - t_begin] = temp_lfp / neuron_list.size();
			}
		}
		GE_file.close();
	} else if (LFP_type == "inh") {
		GI_file.open(dir + "/GI.bin", ios::binary);
		GI_file.read((char*)buffer_vec.data(), 2*sizeof(size_t));

		// For t = [0, t_begin];
		V_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
		GI_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
		// For t = (t_begin, t_end]
		double dym_val[3] = {0.0, 0.0, 0.0};
		double temp_lfp;
		for (size_t i = t_begin; i < t_end; i++) {
			temp_lfp = 0;
			if (neuron_list.size() == shape[1]) {
				for (int j = 0; j < neuron_list.size(); j ++) {
					V_file.read((char*)&dym_val[0], sizeof(double));
					GI_file.read((char*)&dym_val[2], sizeof(double));
					temp_lfp += (- dym_val[gi_id] * (dym_val[v_id] - V_i)) * spatial_weights[j];
				}
				lfp[i - t_begin] = temp_lfp / neuron_list.size();
			} else {
				for (int j = 0; j < diff_list.size() - 1; j ++) {
					V_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
					GI_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
					V_file.read((char*)&dym_val[0], sizeof(double));
					GI_file.read((char*)&dym_val[2], sizeof(double));
					temp_lfp += (- dym_val[gi_id] * (dym_val[v_id] - V_i)) * spatial_weights[j];
				}
				V_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
				GI_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
				lfp[i - t_begin] = temp_lfp / neuron_list.size();
			}
		}
		GI_file.close();
	} else throw runtime_error("ERROR: wrong LFP type (Accessible types: tot, lek, exi, inh.)");
	V_file.close();
}
