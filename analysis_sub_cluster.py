import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import subprocess
import sys
import mylib
import itertools
import random

NN = 11
target_index = 55
loading_dir = '/media/kyle/Drive/ResearchData/Apr22/t3/'
output_dir = './results/Apr22/t3mix-' + str(target_index) + '/'
os.system('mkdir -p ' + output_dir)

# all neurons in the pool
conMat = mylib.load_matrix(loading_dir = loading_dir, filename = 'conMat.txt')
connecting_pool = conMat[target_index,:].nonzero()
connecting_pool = connecting_pool[0].astype(int)
num_in_pool = len(connecting_pool)
total_list = range(100)
nonconnecting_pool = np.array(list(set(total_list) - set(connecting_pool)))
nonconnecting_pool = random.sample(nonconnecting_pool, 5)
mix_pool = np.array(list(connecting_pool) + list(nonconnecting_pool))
num_in_pool = len(mix_pool)

# preparing storage for data;
data_dic = {'index':[], 'connections':[], 'connecting portion':[], 'mean firing rate':[], 'signal noise ratio':[], 'peak time':[], 'time constant':[]}
data_out = pd.DataFrame(data_dic)
counter = 0
for i in range(NN):
    num_in_lfp = i + 1
    if num_in_pool < num_in_lfp:
    	print 'Target number of neurons exceed the total number of neuron in the connected pool.'
        data_out.to_csv(output_dir + "asc_temp.csv", sep = '\t', float_format = '%.6f', index  = False, columns = ['index', 'connections', 'connecting portion', 'mean firing rate', 'signal noise ratio', 'peak time', 'time constant'])
    	sys.exit(1)
    # randomly choose 'number_in_lfp' neurons in the pool and list all possible combinaitions;
    trials = list(itertools.combinations(mix_pool, num_in_lfp))
    num_of_trials = len(trials)

    trials_save = pd.DataFrame(trials)
    trials_save.to_csv(output_dir + 'trials'+ str(num_in_lfp) +'.csv', header = False, index = False)

    # loading neuronal type
    neurons = pd.read_csv(loading_dir + 'postNeuron.txt', sep = '\t', header = None)
    types = neurons[:][0].astype(int)

    time_lb = 1000
    time_ub = 10000
    negative_time_delay = 60
    positive_time_delay = 100
    total_neuron_number = 100
    dt = 0.25
    for j in range(num_of_trials):
		# print '(' + str(i) + ',' +str(j) + ')'
        data_out.loc[counter, 'index'] = counter
        data_out.loc[counter, 'connections'] = num_in_lfp
        # num_exc, num_inh = mylib.DivideNeuronalTypes(neuron_types = types, neuron_list = trials[j])
        portion_counter = 0
        for k in range(num_in_lfp):
            if trials[j][k] in connecting_pool:
                portion_counter += 1
        data_out.loc[counter, 'connecting portion'] = portion_counter * 1.0 / num_in_lfp
        rate = 0
        for k in range(num_in_lfp):
			rate += mylib.mean_rate(loading_dir = loading_dir, filename = 'rasterPost.txt', index = trials[j][k], tmax = 10)
        data_out.loc[counter, 'mean firing rate'] = rate / num_in_lfp
        indices = trials[j]
        indices_str = [str(nn) for nn in indices]
        indices_str = ','.join(indices_str)
        subprocess.call(['./lfp/calculate-lfp.out', loading_dir, str(target_index), indices_str, str(time_lb), str(time_ub), str(total_neuron_number)])
        subprocess.call(['./tdmi/calculate-tdmi.out', str(dt), str(negative_time_delay), str(positive_time_delay)])
        time_series, signal_order, signal_rand = mylib.import_tdmi()
        data_out.loc[counter, 'signal noise ratio'], data_out.loc[counter, 'peak time'], data_out.loc[counter, 'time constant'] = mylib.tdmi_parameters(time_series, signal_order, signal_rand)
        # Draw figures:
        # saving_filename = 'tdmi-' + str(num_in_lfp) + '-' + str(j) + '-' + str(time_lb) + '_' + str(time_ub) + '-';
        # str_dt = str(dt).split('.')
        # saving_filename += str_dt[0] + '_' + str_dt[1] + '-' + str(negative_time_delay) + '-' + str(positive_time_delay)
		# # saving_filename = 'tdmi-20-' + str(num)
		# # print figure_text
        # mylib.PlotTdmi(time_series, signal_order, signal_rand, saving_dir = output_dir, saving_filename = saving_filename)
        counter += 1
        print '================================================'
data_out.to_csv(output_dir + "asc.csv", sep = '\t', float_format = '%.6f', index  = False, columns = ['index', 'connections', 'connecting portion', 'mean firing rate', 'signal noise ratio', 'peak time', 'time constant'])
