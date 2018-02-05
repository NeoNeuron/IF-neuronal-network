//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-08 12:34:37
//	Description: Source code of connectivity_matrix.cpp;
//***************
#include "../include/connectivity_matrix.h"
#include "../include/io.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <numeric>
using namespace std;

void ConnectivityMatrix::Scan(int target_value, int row_index, vector<int> &output_indices) {
	for (int s = 0; s < neuron_number_; s++) {
		if (matrix_[row_index][s] == target_value) output_indices.push_back(s);
	}
}

void ConnectivityMatrix::FindLeastPath() {
	// initialize path_matrix_
	for (int i = 0; i < neuron_number_; i++) {
		for (int j = 0; j < neuron_number_; j++) {
			if (i == j) continue;
			else {
				if (matrix_[i][j] == 1) path_matrix_[i][j] = 1;
				else path_matrix_[i][j] = 1000000;
			}
			mediate_mode_matrix_[i][j] = -1;
		}
	}
	// Floid algorithm for least path length
	for (int k = 0; k < neuron_number_; k++) {
		for (int i = 0; i < neuron_number_; i++) {
			if (i != k) {
				for (int j = 0; j < neuron_number_; j++) {
					if (j!=k && j!=i && path_matrix_[i][j]>path_matrix_[i][k] + path_matrix_[k][j]) {
						path_matrix_[i][j] = path_matrix_[i][k] + path_matrix_[k][j];
						mediate_mode_matrix_[i][j] = k;
					}
				}
			}
		}
	}
}

void ConnectivityMatrix::CalculateClusteringCoefficient() {
	vector<int> neighbor;
	double counter = 0;
	for (int k = 0; k < neuron_number_; k++) {
		for (int s = 0; s < neuron_number_; s++) {
			if (k == s || matrix_[k][s] == 1) neighbor.push_back(s);
		}
		for (int i = 0; i < neighbor.size(); i++) {
			for (int j = 0; j < neighbor.size(); j++) {
				counter += matrix_[neighbor[i]][neighbor[j]];
			}
		}
		clustering_coefficient_[k] = counter / (neighbor.size()*(neighbor.size() - 1));
		neighbor.clear();
		counter = 0;
	}
}

void ConnectivityMatrix::SetNeuronNumber(int neuron_number) {
	neuron_number_ = neuron_number;
	vector<int> add_int0(neuron_number_, 0);
	vector<int> add_intn1(neuron_number_, -1);
	matrix_.clear();
	path_matrix_.clear();
	mediate_mode_matrix_.clear();
	mediate_mode_matrix_.clear();
	clustering_coefficient_.clear();
	matrix_.resize(neuron_number_, add_int0);
	path_matrix_.resize(neuron_number_, add_int0);
	mediate_mode_matrix_.resize(neuron_number_, add_intn1);
	clustering_coefficient_.resize(neuron_number_, 0);
}

void ConnectivityMatrix::SetConnectingDensity(int connecting_density) {
	connecting_density_ = connecting_density; // number of neurons that directly connected with it on its one side;
	for (int i = 0; i < neuron_number_; i++)  {
		for (int j = 0; j < neuron_number_; j++) {
			if (i != j) {
				if (abs(i - j) <= connecting_density_ or neuron_number_ - abs(i - j) <= connecting_density_) {
				matrix_[i][j] = 1;
				path_matrix_[i][j] = 1;
				} else path_matrix_[i][j] = 1000000;
			}
		}
	}
	FindLeastPath();
	CalculateClusteringCoefficient();
}


void ConnectivityMatrix::LoadMatrix(vector<vector<int> >& matrix) {
	for (int i = 0; i < neuron_number_; i++) {
		for (int j = 0; j < neuron_number_; j++) {
			if (i == j) { continue; }
			else {
				matrix_[i][j] = matrix[i][j];
				if (matrix_[i][j] == 1) {
					path_matrix_[i][j] = 1;
				} else {
					path_matrix_[i][j] = 1000000;
				}
			}
		}
	}
}

void ConnectivityMatrix::Rewire(double p, int seed, bool OutputOption) {
	cout << 2 * neuron_number_ * connecting_density_ << " connections total with ";
	srand(seed);
	double x; // random variable;
	int ind, empty_connection, count = 0;
	vector<int> ones, zeros;
	for (int i = 0; i < neuron_number_; i++) {
		Scan(1, i, ones);
		for (int j = 0; j < ones.size(); j++) {
			x = rand() / (RAND_MAX + 1.0);
			if (x <= p) {
				Scan(0, i, zeros);
				for (vector<int>::iterator it = zeros.begin(); it != zeros.end(); it++) {
					if (*it == i) {
						zeros.erase(it);
						break;
					}
				}
				empty_connection = zeros.size();
				ind = rand() % empty_connection;
				matrix_[i][zeros[ind]] = 1;
				matrix_[i][ones[j]] = 0;
				zeros.clear();
				count += 1;
			}
		}
		ones.clear();
	}
	cout << count << " rewirings." << endl;
	if (OutputOption == true) {
		FindLeastPath();
		CalculateClusteringCoefficient();
	}
}

int ConnectivityMatrix::ReadMatrix(int i, int j) {
	return matrix_[i][j];
}

double ConnectivityMatrix::GetMeanPath() {
	double counter = 0;
	for (int i = 0; i < neuron_number_; i++) {
		for (int j = 0; j < neuron_number_; j++) {
			counter += path_matrix_[i][j];
		}
	}
	return counter / (neuron_number_*(neuron_number_ - 1));
}

double ConnectivityMatrix::GetMeanClusteringCoefficient() {
	double counter = 0;
	for (int i = 0; i < neuron_number_; i++) {
		counter += clustering_coefficient_[i];
	}
	return counter / neuron_number_;
}

void ConnectivityMatrix::OutMatrix(string path) {
	Print2DBin(path, matrix_, "trunc");
}

bool ConnectivityMatrix::IsConnect(){
	int counter = 0;
	for (vector<vector<int> >::iterator it = matrix_.begin(); it != matrix_.end(); it ++) {
		counter += accumulate(it->begin(), it->end(), 0);
	}
	if (counter == 0) return false;
	else return true;
}
