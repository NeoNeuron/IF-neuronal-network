//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-02-21 16:05:24
//	Description: source code of multi_network.cpp;
//***************
#include "multi_network.h"
#include <cmath>

using namespace std;

void UpdateSystemState(NeuronalNetwork & pre_network, NeuronalNetwork & post_network, vector<vector<bool> > & connectivity_matrix, double t, double dt) {
	//	Update pre_network in two network system;
	pre_network.UpdateNetworkState(t, dt);
	//	Transmit spiking information from pre_network to post_network;
	vector<vector<Spike> > tempPreSpikes, tempPostSpikes;
	pre_network.OutNewSpikes(t, tempPreSpikes);
	tempPostSpikes.resize(connectivity_matrix[1].size());
	for (int i = 0; i < connectivity_matrix.size(); i++) {	
		if (tempPreSpikes[i].size() != 0) {
			for (int j = 0; j < connectivity_matrix.size(); j++) {
				if (connectivity_matrix[i][j] == true) {
					tempPostSpikes[j].insert(tempPostSpikes[j].end(), tempPreSpikes[i].begin(), tempPreSpikes[i].end());
				}
			}
		}
	}	
	post_network.InNewSpikes(tempPostSpikes);
	//	Update post_network in two network system;
	post_network.UpdateNetworkState(t, dt);
}