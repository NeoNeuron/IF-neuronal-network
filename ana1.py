import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
import os
import subprocess
import mylib

num_in_lfp = 1
loading_dir = '/media/kyle/Drive/ResearchData/Apr22/t3/'

target_index = 80

conMat = mylib.load_matrix(loading_dir = loading_dir, filename = 'conMat.txt')

connecting_pool = conMat[target_index,:].nonzero()

connecting_pool = connecting_pool[0]

num_in_pool = len(connecting_pool)

# loading neuronal type
neurons = pd.read_csv(loading_dir + 'postNeuron.txt', sep = '\t', header = None)
types = neurons[:][0]
# print types[0]

time_lb = 1000
time_ub = 10000
expected_occupancy = 50
negative_time_delay = 60
positive_time_delay = 100
total_neuron_number = 100
num = 0
dt = 0.25
# preparing storage for data;
data_dic = {'index':np.zeros(num_in_pool),'type':np.zeros(num_in_pool), 'mean firing rate':np.zeros(num_in_pool), 'signal noise ratio':np.zeros(num_in_pool), 'peak time':np.zeros(num_in_pool), 'time constant':np.zeros(num_in_pool)}
data_out = pd.DataFrame(data_dic) 	

for i in range(num_in_pool):
	data_out.ix[i]['index'] = int(connecting_pool[i])
	data_out.ix[i]['type'] = int(types[i])
	data_out.ix[i]['mean firing rate'] = mylib.mean_rate(loading_dir = loading_dir, filename = 'rasterPost.txt', index = connecting_pool[i], tmax = 10)
	subprocess.call(["./lfp/calculate-lfp.out", loading_dir, str(int(target_index)), '1', str(time_lb), str(time_ub), 'all', str(i), str(total_neuron_number)])
	subprocess.call(['./tdmi/calculate-tdmi.out', str(expected_occupancy), str(dt), str(negative_time_delay), str(positive_time_delay)])
	data_out.ix[i]['signal noise ratio'], data_out.ix[i]['peak time'], data_out.ix[i]['time constant'] = mylib.tdmi_parameters()

	saving_filename = 'tdmi-' + str(int(connecting_pool[i])) + '-1-' + str(time_lb) + '_' + str(time_ub) + '-all-' + str(expected_occupancy) + '-';
	str_dt = str(dt).split('.')
	saving_filename += str_dt[0] + '_' + str_dt[1] + '-' + str(negative_time_delay) + '-' + str(positive_time_delay)
	# saving_filename = 'tdmi-20-' + str(num)
	# create figure_text
	figure_text = mylib.CreateText(loading_dir = loading_dir, neuron_index = target_index, order = 1, classification = 'all', num = int(connecting_pool[i]))
	# print figure_text
	mylib.PlotTdmi(saving_filename = saving_filename, figure_text = figure_text)

data_out.to_csv(loading_dir + "ana1.csv", float_format = '%.6f', index  = False)