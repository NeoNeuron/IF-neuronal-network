import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
import os
import mylib as my

num_in_lfp = 1
loading_dir = '/media/kyle/Drive/ResearchData/Apr22/t3/'

target_index = 0

conMat = my.load_matrix(loading_dir = loading_dir, filename = 'conMat.txt')

connecting_pool = conMat[target_index,:].nonzero()

connecting_pool = connecting_pool[0]

num_in_pool = len(connecting_pool)

# loading neuronal type
neurons = pd.read_csv(loading_dir + 'preNeuron.txt', sep = '\t', header = None)
types = neurons[:][0]
# print types[0]

time_lb = 1000
time_ub = 10000
expected_occupancy = 50
negative_time_delay = 60
positive_time_delay = 100
num = 0
timing_step_list = 0.25
# preparing storage for data;
data_dic = {'index':np.zeros(num_in_pool),'type':np.zeros(num_in_pool), 'mean firing rate':np.zeros(num_in_pool), 'signal noise ratio':np.zeros(num_in_pool), 'peak time':np.zeros(num_in_pool), 'decay constant':np.zeros(num_in_pool)}
data_out = pd.DataFrame(data_dic) 	

for i in range(num_in_pool):
	data_out.ix[i]['index'] = connecting_pool[i]
	data_out.ix[i]['type'] = types[i]
	