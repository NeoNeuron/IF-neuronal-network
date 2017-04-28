import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
import os
import subprocess
import sys
import mylib
import itertools

num_in_lfp = 8
loading_dir = '/media/kyle/Drive/ResearchData/Apr22/t3/'

target_index = 0
# all neurons in the pool
conMat = mylib.load_matrix(loading_dir = loading_dir, filename = 'conMat.txt')
connecting_pool = conMat[target_index,:].nonzero()
connecting_pool = connecting_pool[0].astype(int)
num_in_pool = len(connecting_pool)
if num_in_pool < num_in_lfp:
	print 'Target number of neurons exceed the total number of neuron in the connected pool.'
	sys.exit(1)
# randomly choose 'number_in_lfp' neurons in the pool and list all possible combinaitions;
trials = list(itertools.combinations(connecting_pool, num_in_lfp))
num_of_trials = len(trials)

trials_save = pd.DataFrame(trials)
trials_save.to_csv('./tdmi/figure-eps/trials.csv', header = False, index = False)

# loading neuronal type
neurons = pd.read_csv(loading_dir + 'postNeuron.txt', sep = '\t', header = None)
types = neurons[:][0].astype(int)
# print types[0]

time_lb = 1000
time_ub = 10000
expected_occupancy = 50
negative_time_delay = 60
positive_time_delay = 100
total_neuron_number = 100
dt = 0.25
# preparing storage for data;
data_dic = {'mean firing rate':np.zeros(num_in_pool), 'signal noise ratio':np.zeros(num_in_pool), 'peak time':np.zeros(num_in_pool), 'time constant':np.zeros(num_in_pool)}
data_out = pd.DataFrame(data_dic) 	

for i in range(num_of_trials):
	rate = 0
	for j in range(num_in_lfp):
		rate += mylib.mean_rate(loading_dir = loading_dir, filename = 'rasterPost.txt', index = trials[i][j], tmax = 10)
	data_out.loc[i, 'mean firing rate'] = rate / num_in_lfp
	indices = trials[i]
	indices_str = [str(nn) for nn in indices]
	indices_str = ','.join(indices_str)
	subprocess.call(['./lfp/calculate-lfp.out', loading_dir, str(target_index), indices_str, str(time_lb), str(time_ub), str(total_neuron_number)])
	subprocess.call(['./tdmi/calculate-tdmi.out', str(expected_occupancy), str(dt), str(negative_time_delay), str(positive_time_delay)])

	time_series, signal_order, signal_rand = mylib.import_tdmi()

	data_out.loc[i, 'signal noise ratio'], data_out.loc[i, 'peak time'], data_out.loc[i, 'time constant'] = mylib.tdmi_parameters(time_series, signal_order, signal_rand)

	saving_filename = 'tdmi-' + str(num_in_lfp) + '-' + str(i) + '-' + str(time_lb) + '_' + str(time_ub) + '-' + str(expected_occupancy) + '-';
	str_dt = str(dt).split('.')
	saving_filename += str_dt[0] + '_' + str_dt[1] + '-' + str(negative_time_delay) + '-' + str(positive_time_delay)
	# saving_filename = 'tdmi-20-' + str(num)
	# print figure_text
	mylib.PlotTdmi(time_series, signal_order, signal_rand, saving_filename = saving_filename)

data_out.to_csv("./tdmi/figure-eps/ana1.csv", sep = '\t', float_format = '%.6f', index  = False, columns = ['mean firing rate', 'signal noise ratio', 'peak time', 'time constant'])