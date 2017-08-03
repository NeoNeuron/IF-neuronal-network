#!/usr/bin/python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import subprocess
import sys

df = float(sys.argv[1]);
duf = float(sys.argv[2]);
frange = sys.argv[3];
frange = [float(i) for i in frange.split(',')]
ufrange = sys.argv[4];
ufrange = [float(i) for i in ufrange.split(',')]

# prepare pathes;
loading_dir = "./tmp/"
saving_dir = './data/Aug03/'

# setting preliminary parameters
total_neuron_number = 100
simulation_accomplish = True
time_lb = 1000
time_ub = 10000
negative_time_delay = 60
positive_time_delay = 100

subprocess.call['./bin/raster.out', './tmp/rasterPre.txt', '0', '1000,20000', 'raster.csv']
subprocess.call['./bin/lfp.out', './tmp/postI.txt', '0', '1000,20000', 'lfp.csv']
subprocess.call['./bin/mi.out', '1', '0.25', '100,200']

subprocess.call['./bin/lfp.out', './tmp/postI.txt', '0,1,2,3,4,5,6,7,8,9', '1000,20000', 'lfp.csv']
subprocess.call['./bin/mi.out', '1', '0.25', '100,200']
subprocess.call['./bin/sta.out', './data/raster/raster.csv', './data/lfp/lfp.csv', '-25,50']
# subprocess.call['./bin/lcc.out', './data/raster/raster.csv', './data/lfp/lfp.csv', '0.25', '100,200']

Neuron cell;
double t = 0, dt = 0.1, tmax = 10000;
double f = f_min, uf = uf_min;
vector<double> impt_e, impt_i;
vector<double> spike_train;
vector<vector<double> > firing_rates;
vector<double> newline;
double u;
while (f <= f_max) {
    newline.clear();
    while (uf <= uf_max) {
        u = uf / f;
        // cout << f << endl;
        cell.SetDrivingType(false);
        cell.SetPoissonRate(true, u);
        cell.SetFeedforwardConductance(true, f);
        // cout << u << ',' << f << endl;
        while (t < tmax) {
            cell.UpdateNeuronalState(t, dt, impt_e, impt_i);
            // cout << v << endl;
            t += dt;
        }
        cell.OutSpikeTrain(spike_train);
        // cout << spike_train.size()*1000.0/tmax << endl;
        newline.push_back(spike_train.size()*1000.0/tmax);
        cell.Reset();
        t = 0;
        uf += duf;
    }
    uf = uf_min;
    f += df;
    firing_rates.push_back(newline);
}

// OUTPUTS:
string path = "./data/tuninng.csv";
Print2D(path, "trunc", firing_rates);


# Generating neruonal data based on settings above;
subprocess.call(["./bin/nets.out", loading_dir])
# Setting loops for local field potentials;
# target_neuron_indice_list = random.sample(all_neuron, 1)
target_neuron_indice_list = range(0, total_neuron_number)
order = 2
# prepare lists:
conMat = mylib.load_matrix(loading_dir = loading_dir, filename = 'conMat.txt')
postMat = mylib.load_matrix(loading_dir = loading_dir, filename = 'postMat.txt')

lists = []
for i in target_neuron_indice_list:
	ll = np.nonzero(conMat[i,:])[0]
	lists.append(ll)
if order == 2:
	lists2 = []
	for i in target_neuron_indice_list:
		ll2 = np.array([]).astype(int)
		# print lists[i]
		for item in lists[i]:
			# print np.nonzero(conMat[item,:])[0]
			ll2 = np.append(ll2, np.nonzero(postMat[item,:])[0])
		# print ll2
		ll2 = np.unique(ll2)
		# print ll2
		# print set(ll2)&set(lists[i])
		ll2 = list(set(ll2) - (set(ll2)&set(lists[i])))
		# print ll2
		lists2.append(ll2)
	lists = lists2

# classification_options = ['all']
# Setting loops for time-delayed mutual information;
timing_step_list = [0.25]
# preparing storage for data;
data_dic = {'index':np.zeros(0).astype(int),'type':np.zeros(0).astype(int), 'mean firing rate':[], 'number of connection':np.zeros(0).astype(int), 'number of excitatory connection':np.zeros(0).astype(int), 'number of inhibitory connection':np.zeros(0).astype(int), 'signal noise ratio':[], 'peak time':[], 'time constant':[]}
data_out = pd.DataFrame(data_dic)

pre_net_types = np.genfromtxt(loading_dir + 'preNeuron.txt', dtype = int, usecols = 0, delimiter = ',')
post_net_types = np.genfromtxt(loading_dir + 'postNeuron.txt', dtype = int, usecols = 0, delimiter = ',')

num_of_trials = 1

# Start loops
for i in range(num_of_trials):
	ind = target_neuron_indice_list[i]
	ll = lists[i]
	# transform list to str
	list_str = [str(nn) for nn in ll]
	list_str = ','.join(list_str)
	# for num in num_list:
	subprocess.call(["./bin/lfp.out", loading_dir, str(int(ind)), list_str, str(time_lb), str(time_ub), str(total_neuron_number)])
	for dt in timing_step_list:
		subprocess.call(['./bin/mi.out', str(dt), str(negative_time_delay), str(positive_time_delay)])
		# create a saving filename
		saving_filename = 'tdmi-' + str(int(ind)) + '-' + str(order) + '-' + str(time_lb) + '-' + str(time_ub) + '-' + str(dt).replace('.', '') + '-' + str(negative_time_delay) + '-' + str(positive_time_delay)

		# create figure_text
		figure_text = mylib.CreateText(loading_dir = loading_dir, neuron_index = int(ind), connecting_list = ll)

		mi_data = pd.read_csv('./data/mi/mi.csv')
		# print figure_text
		mylib.PlotTdmi(time_series = mi_data['timelag'], signal_order = mi_data['ordered'], signal_rand = mi_data['random'], saving_dir = saving_dir, saving_filename = saving_filename, figure_text = figure_text)
		# Output data:
		data_out.loc[ind, 'index'] = ind
		data_out.loc[ind, 'type'] = pre_net_types[ind]
		data_out.loc[ind, 'mean firing rate'] = mylib.mean_rate(loading_dir = loading_dir, filename = 'rasterPre.txt', index = ind, tmax = 10)
		data_out.loc[ind, 'number of connection'] = len(ll)
		data_out.loc[ind, 'number of excitatory connection'], data_out.loc[ind, 'number of inhibitory connection'] = mylib.DivideNeuronalTypes(neuron_types = post_net_types, neuron_list = ll)

		data_out.loc[ind, 'signal noise ratio'],	data_out.loc[ind, 'peak time'], data_out.loc[ind, 'time constant'] = mylib.tdmi_parameters(mi_data)
		print '=================================================='

# print data_2d
data_out.to_csv(saving_dir + "pre-net-data"+str(order)+".csv", float_format = '%.4f', index  = False, columns = ['index','type', 'mean firing rate', 'number of connection', 'number of excitatory connection', 'number of inhibitory connection', 'signal noise ratio', 'peak time', 'time constant'])
