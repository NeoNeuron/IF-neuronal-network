//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-02-21 16:03:06
//	Description: define functions used in multi-network simulation;
//***************
#ifndef _MULTI_NETWORK_H_
#define _MULTI_NETWORK_H_

#include "neuron.h"
#include "connectivity_matrix.h"
#include "group.h"
#include <iostream>
#include <vector>

using namespace std;

void UpdateSystemState(NeuronalNetwork & pre_network, 
												NeuronalNetwork & post_network, 
												vector<vector<bool> > &connectivity_matrix, 
												double t, double dt);

#endif // _MULTI_NETWORK_H_